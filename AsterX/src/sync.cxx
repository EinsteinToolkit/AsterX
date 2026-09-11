#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "../../../CarpetX/CarpetX/src/fillpatch.hxx"
#include "../../../CarpetX/CarpetX/src/schedule.hxx"
#include "../../../CarpetX/CarpetX/src/task_manager.hxx"

namespace AsterX {
using namespace CarpetX;

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

extern "C" void AsterX_RestrictFluxes(CCTK_ARGUMENTS) {
  static const std::vector<int> restrict_groups = {
      CCTK_GroupIndex("AsterX::flux_x"), CCTK_GroupIndex("AsterX::flux_y"),
      CCTK_GroupIndex("AsterX::flux_z")};

  active_levels->loop_fine_to_coarse([&](const auto &leveldata) {
    // Only restrict from a child level that is inside the active window
    // [min_level, max_level), i.e. one that is time-aligned with this level.
    // Under subcycling a coarse-only batch has no active child, so nothing
    // is restricted; without subcycling every level is active, so this is
    // the same as restricting from every level but the finest.
    if (leveldata.level + 1 < active_levels->max_level)
      RestrictNoPoison(cctkGH, leveldata.level, restrict_groups);
  });
}

extern "C" void AsterX_RestrictAuxTermsForAvecPsiRHS(CCTK_ARGUMENTS) {
  static const std::vector<int> restrict_groups = {
      CCTK_GroupIndex("AsterX::G"), CCTK_GroupIndex("AsterX::Ex"),
      CCTK_GroupIndex("AsterX::Ey"), CCTK_GroupIndex("AsterX::Ez")};

  active_levels->loop_fine_to_coarse([&](const auto &leveldata) {
    // Only restrict from a child level that is inside the active window
    // [min_level, max_level), i.e. one that is time-aligned with this level.
    // Under subcycling a coarse-only batch has no active child, so nothing
    // is restricted; without subcycling every level is active, so this is
    // the same as restricting from every level but the finest.
    if (leveldata.level + 1 < active_levels->max_level)
      RestrictNoPoison(cctkGH, leveldata.level, restrict_groups);
  });
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
