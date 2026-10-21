#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuMemory.H>
#include <loop_device.hxx>

#include <cmath>

#include "setup_eos.hxx"
#include "atmo.hxx"

namespace ID_TabEOS_HydroQuantities {
using namespace Loop;

extern "C" void ID_TabEOS_HydroQuantities_TestYe(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_TestYe;
  DECLARE_CCTK_PARAMETERS;

  const auto eos_1p_poly = EOSX::global_eos_1p_poly;
  const auto eos_3p_tab3d = EOSX::global_eos_3p_tab3d;
  amrex::Gpu::DeviceScalar<unsigned int> failures(0);
  auto *failed = failures.dataPtr();

  // Check before final conversion can mask an inconsistent atmosphere Ye.
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL radial_distance =
            sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        const auto atmo = Con2PrimFactory::make_atmo(
            eos_1p_poly, eos_3p_tab3d, radial_distance, rho_abs_min,
            p_atmo, t_atmo, Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo,
            n_temp_atmo, atmo_tol, true, false, Ye_atmo_beq);
        const CCTK_REAL yeL = Ye(p.I);
        if (!std::isfinite(yeL) || yeL < eos_3p_tab3d->rgye.min ||
            yeL > eos_3p_tab3d->rgye.max ||
            (!(rho(p.I) > atmo.rho_cut) && yeL != atmo.ye_atmo))
          amrex::HostDevice::Atomic::Add(failed, 1U);
      });

  amrex::Gpu::streamSynchronize();
  if (failures.dataValue())
    CCTK_ERROR("ID_TabEOS_HydroQuantities test: inconsistent initial Ye");
}

extern "C" void ID_TabEOS_HydroQuantities_TestTempEnt(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_TestTempEnt;
  DECLARE_CCTK_PARAMETERS;

  const auto eos_1p_poly = EOSX::global_eos_1p_poly;
  const auto eos_3p_tab3d = EOSX::global_eos_3p_tab3d;
  const bool temperature_id =
      CCTK_EQUALS(id_temp_ent_type, "constant temperature");
  amrex::Gpu::DeviceScalar<unsigned int> failures(0);
  auto *failed = failures.dataPtr();

  // Check before final conversion can mask an initialization error.
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL radial_distance =
            sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        const auto atmo = Con2PrimFactory::make_atmo(
            eos_1p_poly, eos_3p_tab3d, radial_distance, rho_abs_min,
            p_atmo, t_atmo, Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo,
            n_temp_atmo, atmo_tol, true, false, Ye_atmo_beq);
        const CCTK_REAL rhoL = rho(p.I);
        const CCTK_REAL yeL = Ye(p.I);
        CCTK_REAL temp_expected;
        CCTK_REAL ent_expected;

        if (rhoL > atmo.rho_cut) {
          if (temperature_id) {
            temp_expected = atmo.temp_atmo;
            ent_expected = eos_3p_tab3d->kappa_from_rho_temp_ye(
                rhoL, temp_expected, yeL);
          } else {
            ent_expected = id_entropy;
            temp_expected = eos_3p_tab3d->temp_from_rho_entropy_ye(
                rhoL, ent_expected, yeL);
          }
        } else {
          temp_expected = atmo.temp_atmo;
          ent_expected = atmo.entropy_atmo;
        }

        if (!std::isfinite(temperature(p.I)) ||
            !std::isfinite(entropy(p.I)) ||
            temperature(p.I) != temp_expected ||
            entropy(p.I) != ent_expected)
          amrex::HostDevice::Atomic::Add(failed, 1U);
      });

  amrex::Gpu::streamSynchronize();
  if (failures.dataValue())
    CCTK_ERROR("ID_TabEOS_HydroQuantities test: inconsistent initial "
               "temperature or entropy");
}

extern "C" void ID_TabEOS_HydroQuantities_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_Test;
  DECLARE_CCTK_PARAMETERS;

  const auto eos_1p_poly = EOSX::global_eos_1p_poly;
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

        const CCTK_REAL radial_distance =
            std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        const auto atmo = Con2PrimFactory::make_atmo(
            eos_1p_poly, eos_3p_tab3d, radial_distance, rho_abs_min,
            p_atmo, t_atmo, Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo,
            n_temp_atmo, atmo_tol, true, false, Ye_atmo_beq);

        // Below rho_max, an output inside the cutoff must be a reset.
        // At rho_max, a ceiling clamp can also enter the cutoff; check only
        // EOS consistency there because the original density is not stored.
        if (rhoL < eos_3p_tab3d->rgrho.max && rhoL <= atmo.rho_cut &&
            (rhoL != atmo.rho_atmo || tempL != atmo.temp_atmo ||
             yeL != atmo.ye_atmo || velx(p.I) != 0.0 ||
             vely(p.I) != 0.0 || velz(p.I) != 0.0)) {
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
