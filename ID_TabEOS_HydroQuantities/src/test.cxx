#include <cctk.h>
#include <cctk_Arguments.h>

#include <AMReX.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuMemory.H>
#include <loop_device.hxx>

#include <cmath>

#include "setup_eos.hxx"

namespace ID_TabEOS_HydroQuantities {
using namespace Loop;

extern "C" void ID_TabEOS_HydroQuantities_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_Test;

  const auto eos_3p_tab3d = EOSX::global_eos_3p_tab3d;
  amrex::Gpu::DeviceScalar<unsigned int> failures(0);
  auto *failed = failures.dataPtr();

  // Check the real conversion output, for both initialization modes and
  // graded atmospheres. No grid values are modified by this check.
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL rhoL = rho(p.I);
        const CCTK_REAL tempL = temperature(p.I);
        const CCTK_REAL yeL = Ye(p.I);
        if (!std::isfinite(rhoL) || !std::isfinite(tempL) ||
            !std::isfinite(yeL) ||
            rhoL < eos_3p_tab3d->rgrho.min ||
            rhoL > eos_3p_tab3d->rgrho.max ||
            tempL < eos_3p_tab3d->rgtemp.min ||
            tempL > eos_3p_tab3d->rgtemp.max ||
            yeL < eos_3p_tab3d->rgye.min ||
            yeL > eos_3p_tab3d->rgye.max) {
          amrex::HostDevice::Atomic::Add(failed, 1U);
          return;
        }

        const CCTK_REAL actual[] = {press(p.I), eps(p.I), entropy(p.I)};
        const CCTK_REAL expected[] = {
            eos_3p_tab3d->press_from_rho_temp_ye(rhoL, tempL, yeL),
            eos_3p_tab3d->eps_from_rho_temp_ye(rhoL, tempL, yeL),
            eos_3p_tab3d->entropy_from_rho_temp_ye(rhoL, tempL, yeL)};
        for (int i = 0; i < 3; ++i)
          if (!std::isfinite(actual[i]) || !std::isfinite(expected[i]) ||
              fabs(actual[i] - expected[i]) >
                  2.0e-11 * fmax(fabs(expected[i]), 1.0e-20))
            amrex::HostDevice::Atomic::Add(failed, 1U);
      });

  amrex::Gpu::streamSynchronize();
  if (failures.dataValue())
    CCTK_ERROR("ID_TabEOS_HydroQuantities test: inconsistent initial state");
}

} // namespace ID_TabEOS_HydroQuantities
