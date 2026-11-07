#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <carpetx.hxx>

#include "sync.hxx"

#include <cassert>
#include <vector>

namespace AsterX {

////////////////////////////////////////////////////////////////////////////////
// Level-window helpers (declared in sync.hxx, except has_aligned_child and
// group_indices)

// `level` has a child that is in the active window, i.e. is time-aligned
// with it.
static bool has_aligned_child(const CarpetX::level_t level) {
  return level.index() + 1 < CarpetX::active_window().level_bounds().max;
}

bool full_cascade() { return CarpetX::active_window().level_bounds().min == 0; }

// The flesh's group indices as CarpetX group indices.
static std::vector<CarpetX::gi_t>
group_indices(const std::vector<int> &groups) {
  std::vector<CarpetX::gi_t> gis;
  gis.reserve(groups.size());
  for (const int gi : groups)
    gis.emplace_back(gi);
  return gis;
}

void RestrictFromAlignedChildren(const std::vector<int> &groups) {
  // The window form restricts every level pair (L, L + 1) with both levels in
  // the active window [min_level, max_level), i.e. every level whose child is
  // time-aligned with it (has_aligned_child), finest pair first. Under
  // subcycling a coarse-only batch has no active child, so nothing is
  // restricted; without subcycling every level is active, so this is the same
  // as restricting from every level but the finest. It restricts exactly
  // `groups`, whatever their restrict tags.
  CarpetX::restrict_no_poison(CarpetX::active_window(), group_indices(groups));
}

void SyncGhostsOnly(const std::vector<int> &groups) {
  // CarpetX::sync_ghosts_only copies the same-level ghosts and nothing else;
  // the outer boundary condition is the driver's separate operation, which
  // requires the interior and the ghosts valid.
  const std::vector<CarpetX::gi_t> gis = group_indices(groups);
  CarpetX::sync_ghosts_only(CarpetX::active_window(), gis);
  CarpetX::apply_outer_boundary(CarpetX::active_window(), gis);
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void AsterX_Sync(CCTK_ARGUMENTS) {
  // do nothing
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

  // AsterX_Con2Prim wrote these groups everywhere, but its result at a ghost
  // point is not trusted: the ghosts are copied from the neighbouring
  // interiors (periodic images included) and the outer boundary condition
  // (symmetries included) is applied on top.
  SyncGhostsOnly(groups);
}

// The validity that every variable of `group` has at time level `tl`. A
// write scope that promises exactly this leaves every validity bit as it
// was and poisons nothing.
static CarpetX::valid_t group_validity(const CarpetX::group_t group,
                                       const CarpetX::tl_t tl) {
  const CarpetX::valid_t valid = CarpetX::valid(group, tl, CarpetX::vi_t{0});
  for (int vi = 1; vi < CarpetX::numvars(group); ++vi)
    if (CarpetX::valid(group, tl, CarpetX::vi_t{vi}) != valid)
      CCTK_VERROR("AsterX: the variables of group %s differ in validity at "
                  "time level %d on refinement level %d, so a single write "
                  "scope cannot keep all of them",
                  CCTK_FullGroupName(group.gi.index()), tl.index(),
                  group.level.index());
  return valid;
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
// SyncGhostsOnly would also apply the outer boundary conditions, so the
// exchange is a bare FillBoundary inside a write scope. The scope
// promises the validity each time level already has: FillBoundary only
// overwrites ghost points with copies of valid interior points, so no bit
// changes.
// Single patch only: ghost points on an inter-patch boundary would need
// MultiPatch_Interpolate.
static void FillGhostsFromNeighbours(const cGH *const cctkGH,
                                     const std::vector<int> &groups) {
  assert(CarpetX::num_patches() == 1);
  for (const int gi : groups) {
    const int ntls = CCTK_ActiveTimeLevelsGI(cctkGH, gi);
    assert(ntls > 0);
    const int sync_tl = ntls > 1 ? ntls - 1 : ntls;
    CarpetX::active_window().loop_serially([&](const CarpetX::level_t level) {
      if (!has_aligned_child(level))
        return;
      const CarpetX::group_t group{level, CarpetX::gi_t{gi}};
      const amrex::Periodicity periodicity =
          CarpetX::geometry(level).periodicity();
      for (int itl = 0; itl < sync_tl; ++itl) {
        const CarpetX::tl_t tl{itl};
        CarpetX::write_t scope =
            CarpetX::write(group, tl, group_validity(group, tl),
                           "AsterX::FillGhostsFromNeighbours");
        amrex::MultiFab &mfab = scope.mfab();
        mfab.FillBoundary(0, mfab.nComp(), mfab.nGrowVect(), periodicity);
      }
    });
  }
  CarpetX::synchronize();
}

extern "C" void AsterX_RestrictFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  static const std::vector<int> restrict_groups = {
      CCTK_GroupIndex("AsterX::flux_x"), CCTK_GroupIndex("AsterX::flux_y"),
      CCTK_GroupIndex("AsterX::flux_z")};

  RestrictFromAlignedChildren(restrict_groups);
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

  RestrictFromAlignedChildren(restrict_groups);
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

  CarpetX::prolongate_ghosts(CarpetX::active_window(), group_indices(groups));
}

extern "C" void AsterX_CommdBstag(CCTK_ARGUMENTS) {
  static const std::vector<int> groups = {CCTK_GroupIndex("AsterX::dBx_stag"),
                                          CCTK_GroupIndex("AsterX::dBy_stag"),
                                          CCTK_GroupIndex("AsterX::dBz_stag")};

  SyncGhostsOnly(groups);
}

extern "C" void AsterX_CommdB(CCTK_ARGUMENTS) {
  static const std::vector<int> groups = {CCTK_GroupIndex("AsterX::dB")};

  SyncGhostsOnly(groups);
}

} // namespace AsterX
