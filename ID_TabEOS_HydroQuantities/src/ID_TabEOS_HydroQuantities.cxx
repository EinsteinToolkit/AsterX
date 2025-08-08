#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>

#include <AMReX.H>

#include <loop_device.hxx>

#include "ID_TabEOS_HydroQuantities.hxx"

#define SQ(X) ((X) * (X))

namespace ID_TabEOS_HydroQuantities {

using namespace amrex;
using namespace EOSX;
using namespace Loop;

enum class TS_ID_t { Temperature, Entropy };

extern "C" void ID_TabEOS_HydroQuantities_initial_Y_e(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_initial_Y_e;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING, "Y_e initialization is ENABLED!");

  auto eos_3p_tab3d = global_eos_3p_tab3d;

  // Open the Y_e file, which should countain Y_e(rho) for the EOS table slice
  FILE *Y_e_file = fopen(Y_e_filename, "r");

  // Check if everything is OK with the file
  if ((Y_e_file = fopen(Y_e_filename, "r")) == NULL) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "File \"%s\" does not exist. ABORTING", Y_e_filename);
  } else {
    // Set nrho
    const CCTK_INT nrho = eos_3p_tab3d->interptable->num_points[0];

    // Set interpolation stencil size
    const int interp_stencil_size = 5;

    Ye_reader *id_ye_reader =
        (Ye_reader *)The_Managed_Arena()->alloc(sizeof *id_ye_reader);
    id_ye_reader->init(nrho, Y_e_file, eos_3p_tab3d->interptable->x[0],
                       interp_stencil_size); // eg. logrho array

    // Close the file
    fclose(Y_e_file);

    // Set Y_e
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL radial_distance = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
          CCTK_REAL rho_atm =
              (radial_distance > r_atmo)
                  ? (rho_abs_min * pow((r_atmo / radial_distance), n_rho_atmo))
                  : rho_abs_min;
          rho_atm = std::max(eos_3p_tab3d->rgrho.min, rho_atm);

          if (rho(p.I) > rho_atm * (1 + atmo_tol)) {
            // Interpolate Y_e(rho_i) at gridpoint i
            const CCTK_REAL Y_eL =
                id_ye_reader->interpolate_1d_quantity_as_function_of_rho(
                    interp_stencil_size, nrho, rho(p.I));
            // Finally, set the Y_e gridfunction
            Ye(p.I) = MIN(MAX(Y_eL, eos_3p_tab3d->interptable->xmin<2>()),
                          eos_3p_tab3d->interptable->xmax<2>());
          } else {
            Ye(p.I) = Ye_atmo;
          }
        });

    The_Managed_Arena()->free(id_ye_reader);
  }
}

// Set initial temperature to be constant everywhere (TODO: add other options)
extern "C" void ID_TabEOS_HydroQuantities_initial_temp_ent(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_initial_temp_ent;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING,
             "Temperature and entropy initialization is ENABLED!");

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
        CCTK_REAL radial_distance =
            std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        CCTK_REAL temp_atm =
            (radial_distance > r_atmo)
                ? (t_atmo * std::pow(r_atmo / radial_distance, n_temp_atmo))
                : t_atmo;
        temp_atm = std::max(eos_3p_tab3d->rgtemp.min, temp_atm);
        CCTK_REAL rho_atm =
            (radial_distance > r_atmo)
                ? (rho_abs_min *
                   std::pow((r_atmo / radial_distance), n_rho_atmo))
                : rho_abs_min;
        rho_atm = std::max(eos_3p_tab3d->rgrho.min, rho_atm);

        CCTK_REAL rhoL = rho(p.I);
        CCTK_REAL yeL = Ye(p.I);

        switch (ts_ID) {
        case TS_ID_t::Temperature: {
          temperature(p.I) = temp_atm;
          CCTK_REAL ent_val =
              eos_3p_tab3d->entropy_from_valid_rho_temp_ye(rhoL, temp_atm, yeL);
          entropy(p.I) = ent_val;
          break;
        }
        case TS_ID_t::Entropy: {
          const CCTK_REAL rho_atmo_cut = rho_atm * (1 + atmo_tol);
          if (rhoL > rho_atmo_cut) {
            CCTK_REAL ent_val = id_entropy;
            CCTK_REAL temp_val = eos_3p_tab3d->temp_from_valid_rho_entropy_ye(
                rhoL, ent_val, yeL);
            entropy(p.I) = ent_val;
            temperature(p.I) = temp_val;
          } else {
            temperature(p.I) = temp_atm;
            CCTK_REAL ent_val = eos_3p_tab3d->entropy_from_valid_rho_temp_ye(
                rhoL, temp_atm, yeL);
            entropy(p.I) = ent_val;
          }
          break;
        }
        default:
          assert(0);
        }

        if (temperature(p.I) < 0.0) {
          printf("Negative input for temperature at I=%d (x=%.5e y=%.5e "
                 "z=%.5e): temp=%.5e\n",
                 p.I, p.x, p.y, p.z, temperature(p.I));
        }
      });
}

// Now recompute all HydroQuantities, to ensure consistent initial data
extern "C" void
ID_TabEOS_HydroQuantities_recompute_HydroBase_variables(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_recompute_HydroBase_variables;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING, "Recomputing all HydroBase quantities ...");

  auto eos_3p_tab3d = global_eos_3p_tab3d;

  // table minimum values
  const CCTK_REAL Tmin = eos_3p_tab3d->rgtemp.min;
  const CCTK_REAL rho_min = eos_3p_tab3d->rgrho.min;
  const CCTK_REAL Ye_min = eos_3p_tab3d->rgye.min;
  const CCTK_REAL Ye_max = eos_3p_tab3d->rgye.max;
  const CCTK_REAL eps_min = eos_3p_tab3d->rgeps.min;

  // compute P_min from table at (rho_min, Tmin, Ye_min)
  const CCTK_REAL P_min =
      eos_3p_tab3d->press_from_valid_rho_temp_ye(rho_min, Tmin, Ye_min);

  // Loop over the grid, recomputing the HydroBase quantities
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Find Atmospheric Density
        CCTK_REAL rhoL = rho(p.I);
        CCTK_REAL tempL = temperature(p.I);
        if (!std::isfinite(tempL) || tempL < Tmin) {
          tempL = Tmin;
          temperature(p.I) = tempL;
        }

        CCTK_REAL radial_distance =
            std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        CCTK_REAL rho_atm =
            (radial_distance > r_atmo)
                ? (rho_abs_min * pow((r_atmo / radial_distance), n_rho_atmo))
                : rho_abs_min;
        rho_atm = std::max(rho_atm, rho_min);
        const CCTK_REAL rho_atmo_cut = rho_atm * (1 + atmo_tol);

        if (rhoL > rho_atmo_cut) {
          CCTK_REAL yeL = Ye(p.I);
          yeL = std::clamp(yeL, Ye_min, Ye_max);
          Ye(p.I) = yeL;

          CCTK_REAL Pval =
              eos_3p_tab3d->press_from_valid_rho_temp_ye(rhoL, tempL, yeL);

          if (!std::isfinite(Pval) || Pval < P_min) {
            Pval = P_min;
          }
          press(p.I) = Pval;

          CCTK_REAL eps_val =
              eos_3p_tab3d->eps_from_valid_rho_temp_ye(rhoL, tempL, yeL);
          if (!std::isfinite(eps_val) || eps_val < eps_min) {
            eps_val = eps_min;
          }
          eps(p.I) = eps_val;

        } else {
          // Reset to atmosphere
          CCTK_REAL temp_atmL = tempL;
          rho(p.I) = rho_atm;
          Ye(p.I) = Ye_atmo;

          CCTK_REAL Pval_atm = eos_3p_tab3d->press_from_valid_rho_temp_ye(
              rho_atm, temp_atmL, Ye_atmo);

          if (!std::isfinite(Pval_atm) || Pval_atm < P_min) {
            Pval_atm = P_min;
          }
          press(p.I) = Pval_atm;

          CCTK_REAL eps_val_atm = eos_3p_tab3d->eps_from_valid_rho_temp_ye(
              rho_atm, temp_atmL, Ye_atmo);

          if (!std::isfinite(eps_val_atm) || eps_val_atm < eps_min) {
            eps_val_atm = eps_min;
          }
          eps(p.I) = eps_val_atm;

          velx(p.I) = 0.0;
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
        }
      });
}

} // namespace ID_TabEOS_HydroQuantities
