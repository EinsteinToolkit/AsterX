#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <algorithm>
#include <cmath>

#include <AMReX.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuMemory.H>

#include <loop_device.hxx>

#include "ID_TabEOS_HydroQuantities.hxx"
#include "atmo.hxx"

#define SQ(X) ((X) * (X))

namespace ID_TabEOS_HydroQuantities {

using namespace amrex;
using namespace Con2PrimFactory;
using namespace EOSX;
using namespace Loop;

enum class TS_ID_t { Temperature, Entropy };

extern "C" void ID_TabEOS_HydroQuantities_initial_Y_e(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_initial_Y_e;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING, "Y_e initialization is ENABLED!");

  auto eos_1p_poly = global_eos_1p_poly;
  auto eos_3p_tab3d = global_eos_3p_tab3d;

  // Open the Y_e file, which should countain Y_e(rho) for the EOS table slice
  FILE *const Y_e_file = fopen(Y_e_filename, "r");

  // Check if everything is OK with the file
  if (Y_e_file == NULL) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "File \"%s\" does not exist. ABORTING", Y_e_filename);
  } else {
    // Set nrho
    const CCTK_INT nrho = eos_3p_tab3d->interptable->num_points[0];

    // Set interpolation stencil size
    const int interp_stencil_size = 5;

    Ye_reader id_ye_reader;
    id_ye_reader.init(nrho, Y_e_file, eos_3p_tab3d->interptable->x[0],
                      interp_stencil_size); // eg. logrho array

    // Close the file
    fclose(Y_e_file);

    // Set Y_e
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const CCTK_REAL radial_distance =
              sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
          const auto atmo = make_atmo(
              eos_1p_poly, eos_3p_tab3d, radial_distance, rho_abs_min,
              p_atmo, t_atmo, Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo,
              n_temp_atmo, atmo_tol, true, false, Ye_atmo_beq);

          if (rho(p.I) > atmo.rho_cut) {
            // Interpolate Y_e(rho_i) at gridpoint i
            const CCTK_REAL Y_eL =
                id_ye_reader.interpolate_1d_quantity_as_function_of_rho(
                    interp_stencil_size, nrho, rho(p.I));
            // Finally, set the Y_e gridfunction
            Ye(p.I) = MIN(MAX(Y_eL, eos_3p_tab3d->interptable->xmin<2>()),
                          eos_3p_tab3d->interptable->xmax<2>());
          } else {
            Ye(p.I) = atmo.ye_atmo;
          }
        });

    // loop_all_device is asynchronous. Keep the reader's managed arrays alive
    // until the device has finished using the by-value reader copy.
    Gpu::synchronize();
    id_ye_reader.release();
  }
}

// Set the initial temperature or entropy profile.
extern "C" void ID_TabEOS_HydroQuantities_initial_temp_ent(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_initial_temp_ent;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING,
             "Temperature and entropy initialization is ENABLED!");

  auto eos_1p_poly = global_eos_1p_poly;
  auto eos_3p_tab3d = global_eos_3p_tab3d;

  TS_ID_t ts_ID;

  if (CCTK_EQUALS(id_temp_ent_type, "constant temperature")) {
    ts_ID = TS_ID_t::Temperature;
  } else if (CCTK_EQUALS(id_temp_ent_type, "constant entropy")) {
    ts_ID = TS_ID_t::Entropy;
  } else {
    CCTK_ERROR("Unknown value for parameter \"id_temp_ent_type\"");
  }

  // Loop over the grid, initializing the temperature
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL radial_distance =
            std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        const auto atmo = make_atmo(
            eos_1p_poly, eos_3p_tab3d, radial_distance, rho_abs_min,
            p_atmo, t_atmo, Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo,
            n_temp_atmo, atmo_tol, true, false, Ye_atmo_beq);

        const CCTK_REAL rhoL = rho(p.I);
        const CCTK_REAL yeL = Ye(p.I);

        if (rhoL > atmo.rho_cut) {
          switch (ts_ID) {
          case TS_ID_t::Temperature:
            temperature(p.I) = atmo.temp_atmo;
            entropy(p.I) = eos_3p_tab3d->kappa_from_rho_temp_ye(
                rhoL, atmo.temp_atmo, yeL);
            break;
          case TS_ID_t::Entropy: {
            CCTK_REAL ent_val = id_entropy;
            temperature(p.I) = eos_3p_tab3d->temp_from_rho_entropy_ye(
                rhoL, ent_val, yeL);
            entropy(p.I) = ent_val;
            break;
          }
          default:
            assert(0);
          }
        } else {
          temperature(p.I) = atmo.temp_atmo;
          entropy(p.I) = atmo.entropy_atmo;
        }

        if (temperature(p.I) < 0.0) {
          printf("Negative input for temperature at (x=%.5e y=%.5e "
                 "z=%.5e): temp=%.5e\n",
                 p.x, p.y, p.z, temperature(p.I));
        }
      });
}

// Now recompute all HydroQuantities, to ensure consistent initial data
extern "C" void
ID_TabEOS_HydroQuantities_recompute_HydroBase_variables(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_recompute_HydroBase_variables;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING, "Recomputing all HydroBase quantities ...");

  auto eos_1p_poly = global_eos_1p_poly;
  auto eos_3p_tab3d = global_eos_3p_tab3d;

  // Table validity bounds
  const CCTK_REAL Tmin = eos_3p_tab3d->rgtemp.min;
  const CCTK_REAL Tmax = eos_3p_tab3d->rgtemp.max;
  const CCTK_REAL rho_min = eos_3p_tab3d->rgrho.min;
  const CCTK_REAL rho_max = eos_3p_tab3d->rgrho.max;
  const CCTK_REAL Ye_min = eos_3p_tab3d->rgye.min;
  const CCTK_REAL Ye_max = eos_3p_tab3d->rgye.max;

  // Device loops cannot call CCTK_ERROR. Report failures on the host.
  amrex::Gpu::DeviceScalar<unsigned int> failures(0);
  auto *failed = failures.dataPtr();

  // Loop over the grid, recomputing the HydroBase quantities
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL rhoL = rho(p.I);
        CCTK_REAL tempL = temperature(p.I);
        CCTK_REAL yeL = Ye(p.I);

        const CCTK_REAL radial_distance =
            std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        // Tabulated initial data uses thermal, temperature-primary atmosphere.
        const auto atmo = make_atmo(
            eos_1p_poly, eos_3p_tab3d, radial_distance, rho_abs_min,
            p_atmo, t_atmo, Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo,
            n_temp_atmo, atmo_tol, true, false, Ye_atmo_beq);

        CCTK_REAL Pval;
        CCTK_REAL eps_val;
        CCTK_REAL ent_val;

        // Preserve the existing comparison, including its NaN-density handling.
        if (!(rhoL > atmo.rho_cut)) {
          // Copy the complete atmosphere; do not apply further thermal floors.
          rhoL = atmo.rho_atmo;
          tempL = atmo.temp_atmo;
          yeL = atmo.ye_atmo;
          Pval = atmo.press_atmo;
          eps_val = atmo.eps_atmo;
          ent_val = atmo.entropy_atmo;
          velx(p.I) = 0.0;
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
        } else {
          // Retain the existing fallback for a non-finite input temperature.
          if (!std::isfinite(tempL))
            tempL = Tmin;
          tempL = std::clamp(tempL, Tmin, Tmax);

          if (!std::isfinite(rhoL) || !std::isfinite(yeL)) {
            amrex::HostDevice::Atomic::Add(failed, 1U);
            return;
          }

          rhoL = std::clamp(rhoL, rho_min, rho_max);
          yeL = std::clamp(yeL, Ye_min, Ye_max);

          // rho, T and Ye are authoritative. Do not independently floor
          // pressure or energy after computing them from this state.
          Pval = eos_3p_tab3d->press_from_rho_temp_ye(rhoL, tempL, yeL);
          eps_val = eos_3p_tab3d->eps_from_rho_temp_ye(rhoL, tempL, yeL);
          // Use evolved kappa and reuse the known temperature.
          ent_val = eos_3p_tab3d->kappa_from_rho_temp_ye(rhoL, tempL, yeL);
        }

        if (!std::isfinite(rhoL) || !std::isfinite(tempL) ||
            !std::isfinite(yeL) || !std::isfinite(Pval) ||
            !std::isfinite(eps_val) || !std::isfinite(ent_val)) {
          amrex::HostDevice::Atomic::Add(failed, 1U);
          return;
        }

        rho(p.I) = rhoL;
        temperature(p.I) = tempL;
        Ye(p.I) = yeL;
        press(p.I) = Pval;
        eps(p.I) = eps_val;
        entropy(p.I) = ent_val;
      });

  // Complete the loop before reading or releasing the temporary flag.
  amrex::Gpu::streamSynchronize();
  if (failures.dataValue())
    CCTK_ERROR("Non-finite rho, Ye or EOS output during initial-data conversion");
}

} // namespace ID_TabEOS_HydroQuantities
