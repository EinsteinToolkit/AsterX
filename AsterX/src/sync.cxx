#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "../../../CarpetX/CarpetX/src/fillpatch.hxx"
#include "../../../CarpetX/CarpetX/src/schedule.hxx"
#include "../../../CarpetX/CarpetX/src/task_manager.hxx"

#include "sync.hxx"

#include <cassert>
#include <vector>

namespace AsterX {
using namespace CarpetX;

////////////////////////////////////////////////////////////////////////////////
// Level-window helpers (declared in sync.hxx)

bool has_aligned_child(const int level) {
  assert(active_levels);
  return level + 1 < active_levels->max_level;
}

bool full_cascade() {
  assert(active_levels);
  return active_levels->min_level == 0;
}

void RestrictFromAlignedChildren(const cGH *const cctkGH,
                                 const std::vector<int> &groups) {
  assert(active_levels);
  active_levels->loop_fine_to_coarse([&](const auto &leveldata) {
    // Only restrict from a child level that is inside the active window
    // [min_level, max_level), i.e. one that is time-aligned with this level.
    // Under subcycling a coarse-only batch has no active child, so nothing
    // is restricted; without subcycling every level is active, so this is
    // the same as restricting from every level but the finest.
    if (has_aligned_child(leveldata.level))
      RestrictNoPoison(cctkGH, leveldata.level, groups);
  });
}

void SyncGhostsOnly(const cGH *const cctkGH, const std::vector<int> &groups) {
  SyncGroupsByDirIGhostOnly(cctkGH, groups.size(), groups.data(), nullptr);
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void AsterX_Sync(CCTK_ARGUMENTS) {
  // do nothing
}

void ApplyOuterBC(CCTK_ARGUMENTS, const std::vector<int> &groups) {
  task_manager tasks1;
  task_manager tasks2;

  for (const int gi : groups) {
    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);

      const int ntls = groupdata.mfab.size();
      const int sync_tl = ntls > 1 ? ntls - 1 : ntls;

      // Copy from adjacent boxes on same level and apply boundary conditions
      // Even though this introduces additional communication, it makes the code
      // compatible with symmetries.
      for (int tl = 0; tl < sync_tl; ++tl) {
        tasks1.submit_serially([&tasks2, &leveldata, &groupdata, tl]() {
          FillPatch_Sync(tasks2, groupdata, *groupdata.mfab.at(tl),
                         ghext->patchdata.at(leveldata.patch)
                             .amrcore->Geom(leveldata.level));
        });
      } // for tl
    });
  } // for gi

  tasks1.run_tasks_serially();
  synchronize();
  tasks2.run_tasks_serially();
  synchronize();

  assert(ghext->num_patches() == 1);
}

extern "C" void AsterX_ApplyOuterBCOnPrim(CCTK_ARGUMENTS) {
  static const std::vector<int> groups = {
      CCTK_GroupIndex("HydroBaseX::rho"),
      CCTK_GroupIndex("HydroBaseX::vel"),
      CCTK_GroupIndex("HydroBaseX::eps"),
      CCTK_GroupIndex("HydroBaseX::press"),
      CCTK_GroupIndex("HydroBaseX::Bvec"),
      CCTK_GroupIndex("HydroBaseX::temperature"),
      CCTK_GroupIndex("HydroBaseX::entropy"),
      CCTK_GroupIndex("HydroBaseX::Ye"),
      CCTK_GroupIndex("AsterX::zvec"),
      CCTK_GroupIndex("AsterX::svec"),
      CCTK_GroupIndex("AsterX::dBx_stag"),
      CCTK_GroupIndex("AsterX::dBy_stag"),
      CCTK_GroupIndex("AsterX::dBz_stag")};

  ApplyOuterBC(CCTK_PASS_CTOC, groups);
}

// Same-level ghost exchange of `groups` on every active level that has an
// aligned child, i.e. on the levels RestrictFromAlignedChildren has just
// rewritten; the other levels were not restricted, so their ghost copies still
// match the neighbouring interiors.
// No outer boundary conditions are applied: ghost points that lie outside the
// domain or beyond a refinement boundary keep their locally computed values,
// only the ghost points covered by a neighbouring box are refreshed from its
// interior. In particular the ghost points beyond a reflection symmetry plane
// are not refreshed (that needs the parities of `groups`, which the flux and G
// groups do not declare), so they keep the unrestricted values even where
// their mirror image was restricted.
// Single patch only: ghost points on an inter-patch boundary would need
// MultiPatch_Interpolate.
static void FillGhostsFromNeighbours(const cGH *const cctkGH,
                                     const std::vector<int> &groups) {
  assert(active_levels);
  assert(ghext->num_patches() == 1);
  for (const int gi : groups) {
    active_levels->loop_serially([&](auto &restrict leveldata) {
      if (!has_aligned_child(leveldata.level))
        return;
      auto &restrict groupdata = *leveldata.groupdata.at(gi);
      const int ntls = groupdata.mfab.size();
      const int sync_tl = ntls > 1 ? ntls - 1 : ntls;
      const auto &geom =
          ghext->patchdata.at(leveldata.patch).amrcore->Geom(leveldata.level);
      for (int tl = 0; tl < sync_tl; ++tl) {
        auto &mfab = *groupdata.mfab.at(tl);
        mfab.FillBoundary(0, mfab.nComp(), mfab.nGrowVect(),
                          geom.periodicity());
      }
    });
  }
  synchronize();
}

extern "C" void AsterX_RestrictFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  static const std::vector<int> restrict_groups = {
      CCTK_GroupIndex("AsterX::flux_x"), CCTK_GroupIndex("AsterX::flux_y"),
      CCTK_GroupIndex("AsterX::flux_z")};

  RestrictFromAlignedChildren(cctkGH, restrict_groups);
  // The restriction only rewrites the coarse valid regions. At
  // hydro_correction_order > 2 the flux-difference stencil in AsterX_RHS reads
  // one face beyond each box's interior, so refresh the same-level ghost copies
  // from the (now restricted) neighbouring interiors; otherwise a cell next to
  // an interprocess boundary sees a restricted flux on one side and the stale
  // unrestricted coarse flux on the other. At 2nd order no ghost face is read.
  if (hydro_correction_order > 2)
    FillGhostsFromNeighbours(cctkGH, restrict_groups);
}

// Scheduled only with restrict_fluxes, without subcycling and without
// freeze_evolution; warns only if CarpetX::do_reflux has been switched on
// (it is off by default). The coarse flux the register accumulates at each
// stage is then already the restricted fine flux, so the register holds
// round-off and the end-of-step reflux changes nothing but costs a register
// update per stage.
extern "C" void AsterX_CheckRedundantReflux(CCTK_ARGUMENTS) {
  // CarpetX::do_reflux and CarpetX::max_num_levels are private to the
  // driver; read them by name
  int do_reflux_type;
  const void *const do_reflux_p =
      CCTK_ParameterGet("do_reflux", "CarpetX", &do_reflux_type);
  assert(do_reflux_p);
  assert(do_reflux_type == PARAMETER_BOOLEAN);
  const CCTK_INT do_reflux = *static_cast<const CCTK_INT *>(do_reflux_p);

  int max_num_levels_type;
  const void *const max_num_levels_p =
      CCTK_ParameterGet("max_num_levels", "CarpetX", &max_num_levels_type);
  assert(max_num_levels_p);
  assert(max_num_levels_type == PARAMETER_INT);
  const CCTK_INT max_num_levels =
      *static_cast<const CCTK_INT *>(max_num_levels_p);

  if (!do_reflux || max_num_levels <= 1)
    return;
  // Every process reaches the same verdict; warn once per run
  if (CCTK_MyProc(nullptr) != 0)
    return;
  CCTK_VWARN(CCTK_WARN_ALERT,
             "AsterX::restrict_fluxes = yes and CarpetX::do_reflux = yes "
             "without subcycling: the fluxes are restricted fine-to-coarse at "
             "every RK stage, so the flux register only accumulates round-off "
             "and the reflux of cons_vector has no effect. Leave "
             "CarpetX::do_reflux at its default (no) to skip the register "
             "work.");
}

extern "C" void AsterX_RestrictAuxTermsForAvecPsiRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  // Only G is restricted. The edge electric field E is not: CalcRHSofAvec_impl
  // reads it pointwise, so a restricted E would only change the Avec RHS on
  // the fine-covered coarse edges themselves, and the driver's restriction of
  // Avec overwrites those edges from the fine level anyway (in the AsterX_Sync
  // of every RK stage with CarpetX::restrict_during_sync = yes). The D_i G
  // stencil instead reaches from uncovered coarse edges next to a refinement
  // boundary onto fine-covered vertices, so restricting G does change the
  // evolved coarse state.
  static const std::vector<int> restrict_groups = {
      CCTK_GroupIndex("AsterX::G")};

  RestrictFromAlignedChildren(cctkGH, restrict_groups);
  // Same as in AsterX_RestrictFluxes: at mag_correction_order > 2 the D_i G
  // stencil in CalcRHSofAvec_impl reads one ghost vertex of G, so its ghost
  // copies must match the restricted interiors of the neighbouring boxes. Only
  // the generalized Lorenz gauge has that term; in the algebraic gauge G is
  // not read at all.
  if (mag_correction_order > 2 &&
      CCTK_EQUALS(vector_potential_gauge, "generalized Lorenz"))
    FillGhostsFromNeighbours(cctkGH, restrict_groups);
}

extern "C" void AsterX_ProlongatedBstag(CCTK_ARGUMENTS) {
  static const std::vector<int> groups = {CCTK_GroupIndex("AsterX::dBx_stag"),
                                          CCTK_GroupIndex("AsterX::dBy_stag"),
                                          CCTK_GroupIndex("AsterX::dBz_stag")};

  SyncGroupsByDirIProlongateOnly(cctkGH, groups.size(), groups.data(), nullptr);
}

extern "C" void AsterX_CommdBstag(CCTK_ARGUMENTS) {
  static const std::vector<int> groups = {CCTK_GroupIndex("AsterX::dBx_stag"),
                                          CCTK_GroupIndex("AsterX::dBy_stag"),
                                          CCTK_GroupIndex("AsterX::dBz_stag")};

  SyncGroupsByDirIGhostOnly(cctkGH, groups.size(), groups.data(), nullptr);
}

extern "C" void AsterX_CommdB(CCTK_ARGUMENTS) {
  static const std::vector<int> groups = {CCTK_GroupIndex("AsterX::dB")};

  SyncGroupsByDirIGhostOnly(cctkGH, groups.size(), groups.data(), nullptr);
}

} // namespace AsterX
