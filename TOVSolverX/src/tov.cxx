/* file    tov.cxx
 * author  Frank Loeffler, converted from fortran thorn by Ian Hawk
 *                         modified to be compatible with CarpetX by Johnny
 *                         original thorn from
 * Cactus/arrangements/EinsteinInitialData/TOVSolverX date    2022/07/24 desc
 * Scheduled functions for TOV initial data compatible with CarpetX grid
 */

#include <loop.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>

#include "tov.hxx"
#include "solve_tov_1d.hxx"

// initial setup adapted from Springel+2010

namespace TOVSolverX {
using namespace Loop;

CCTK_REAL *TOV_Surface = 0;
CCTK_REAL *TOV_R_Surface = 0;
CCTK_REAL *TOV_RProp_Surface = 0;

CCTK_REAL *TOV_r_1d = 0;
CCTK_REAL *TOV_rbar_1d = 0;
CCTK_REAL *TOV_press_1d = 0;
CCTK_REAL *TOV_phi_1d = 0;
CCTK_REAL *TOV_m_1d = 0;
CCTK_REAL *TOV_mbary_1d = 0;
CCTK_REAL *TOV_rprop_1d = 0;

/*@@
   @routine    TOVX_Exact
   @date       Fri Jul 29 10:00:00 2022
   @author     Frank Loeffler, converted fortran routine by Ian Hawke
   @desc
       Schedule routine for interpolation of 1D to 3D grid
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
@@*/
extern "C" void TOVX_C_Exact(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TOVX_C_Exact;
  DECLARE_CCTK_PARAMETERS;

  assert(TOV_Surface != 0);
  assert(TOV_R_Surface != 0);
  assert(TOV_r_1d != 0);
  assert(TOV_rbar_1d != 0);
  assert(TOV_press_1d != 0);
  assert(TOV_phi_1d != 0);
  assert(TOV_m_1d != 0);

  /* allocate local arrays */
  CCTK_REAL *r_to_star = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));
  CCTK_REAL *press_point = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));
  CCTK_REAL *rho_point = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));
  CCTK_REAL *eps_point = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));
  CCTK_REAL *mu_point = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));
  CCTK_REAL *phi_point = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));
  CCTK_REAL *r_point = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));

  /* clear initial data */
  if (TOV_Clear_Initial_Data > 0 && !(TOV_Use_Old_Initial_Data)) {
    grid.loop_int<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          alp_cell(p.I) = 1.0;
          betax_cell(p.I) = 0.0;
          betay_cell(p.I) = 0.0;
          betaz_cell(p.I) = 0.0;

          dtalp_cell(p.I) = 0.0;
          dtbetax_cell(p.I) = 0.0;
          dtbetay_cell(p.I) = 0.0;
          dtbetaz_cell(p.I) = 0.0;

          gxx_cell(p.I) = 0.0;
          gyy_cell(p.I) = 0.0;
          gzz_cell(p.I) = 0.0;
          gxy_cell(p.I) = gxz_cell(p.I) = gyz_cell(p.I) = 0.0;

          kxx_cell(p.I) = 0.0;
          kyy_cell(p.I) = 0.0;
          kzz_cell(p.I) = 0.0;
          kxy_cell(p.I) = kxz_cell(p.I) = kyz_cell(p.I) = 0.0;
        });
  }

  if (!TOV_Use_Old_Matter_Initial_Data) {
    grid.loop_int<1, 1, 1>(grid.nghostzones,
                           [=] CCTK_HOST(const Loop::PointDesc &p)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 rho(p.I) = 0.0;
                                 press(p.I) = 0.0;
                                 eps(p.I) = 0.0;
                                 velx(p.I) = 0.0;
                                 vely(p.I) = 0.0;
                                 velz(p.I) = 0.0;
                               });

    grid.loop_int<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_int<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_int<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });
  }

  CCTK_INT tov_lapse = CCTK_EQUALS(initial_lapse, "tov");
  CCTK_INT tov_shift = CCTK_EQUALS(initial_shift, "tov");

  tov_combine_method_t combine_method;
  if (CCTK_EQUALS(TOV_Combine_Method, "maximum")) {
    combine_method = tov_combine_method_t::maximum;
  } else if (CCTK_EQUALS(TOV_Combine_Method, "average")) {
    combine_method = tov_combine_method_t::average;
  } else {
    CCTK_ERROR("Unknown value for parameter \"TOV_Combine_Method\"");
  }

  auto TOV_press_1d_local = TOV_press_1d;
  auto TOV_phi_1d_local = TOV_phi_1d;
  auto TOV_rbar_1d_local = TOV_rbar_1d;
  auto TOV_r_1d_local = TOV_r_1d;
  auto TOV_Surface_local = TOV_Surface;
  auto TOV_R_Surface_local = TOV_R_Surface;

  grid.loop_int<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_INT star;
        CCTK_REAL g_diag, max_g_diag, max_rho;
        CCTK_REAL my_velx, my_vely, my_velz;
        CCTK_REAL local_tiny = 1.0e-14;
        CCTK_REAL my_psi4 = 1.0;

        for (star = 0; star < TOV_Num_TOVs; star++) {
          r_to_star[star] =
              sqrt((p.x - TOV_Position_x[star]) * (p.x - TOV_Position_x[star]) +
                   (p.y - TOV_Position_y[star]) * (p.y - TOV_Position_y[star]) +
                   (p.z - TOV_Position_z[star]) * (p.z - TOV_Position_z[star]));
          int star_i = star * TOV_Num_Radial;

          /* do the actual interpolation */
          TOVX_C_interp_tov_isotropic(
              star, TOV_Num_Radial, TOV_Fast_Interpolation,
              &(TOV_press_1d_local[star_i]), &(TOV_phi_1d_local[star_i]),
              &(TOV_rbar_1d_local[star_i]), &(TOV_r_1d_local[star_i]),
              &(r_to_star[star]), TOV_Surface_local[star], &(press_point[star]),
              &(phi_point[star]), &(r_point[star]));

          /* is some perturbation wanted? */
          if (Perturb[star] == 0)
            rho_point[star] = pow(press_point[star] / TOV_K, 1.0 / TOV_Gamma);
          else {
            CCTK_REAL r = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
            CCTK_REAL PI = 4.0 * atan(1.0);
            rho_point[star] =
                pow(press_point[star] / TOV_K, 1.0 / TOV_Gamma) *
                (1.0 + Pert_Amplitude[star] *
                           cos(PI / 2.0 * r / TOV_R_Surface_local[star]));
          }
          /* set eps_point */
          if (rho_point[star] > local_tiny)
            eps_point[star] =
                press_point[star] / (TOV_Gamma - 1.0) / rho_point[star];
          else
            eps_point[star] = 0.0;
          /* set mu_point */
          mu_point[star] = rho_point[star] * (1.0 + eps_point[star]);
        }

        /* find out from which star we want to have the data */
        if (combine_method == tov_combine_method_t::maximum) {
          star = 0;
          max_g_diag = 0.0;
          max_rho = rho_point[0];
          for (int star_i = 0; star_i < TOV_Num_TOVs; star_i++) {
            g_diag = (r_point[star_i] / (r_to_star[star_i] + 1.0e-30)) *
                     (r_point[star_i] / (r_to_star[star_i] + 1.0e-30));
            if ((g_diag - max_g_diag) > local_tiny) {
              max_g_diag = g_diag;
              star = star_i;
            }
            if ((rho_point[star_i] - max_rho) > local_tiny) {
              max_rho = rho_point[star_i];
              star = star_i;
            }
          }

          /* handle initial data */
          if (TOV_Use_Old_Initial_Data) {
            /* check metric */
            if ((my_psi4 * gxx_cell(p.I) < max_g_diag) &&
                (my_psi4 * gyy_cell(p.I) < max_g_diag) &&
                (my_psi4 * gzz_cell(p.I) < max_g_diag)) {
              gxx_cell(p.I) = max_g_diag / my_psi4;
              gyy_cell(p.I) = max_g_diag / my_psi4;
              gzz_cell(p.I) = max_g_diag / my_psi4;
              gxy_cell(p.I) = gxz_cell(p.I) = gyz_cell(p.I) = 0.0;
            }
            /* check matter */
            if (TOV_Use_Old_Matter_Initial_Data) {
              if (rho(p.I) > max_rho) {
                /* we do not need this array element anymore, since we use
                 * the initial data, so lets use it */
                star = 0;
                max_rho = rho(p.I);
                eps_point[star] = eps(p.I);
                press_point[star] = press(p.I);
                my_velx = velx(p.I);
                my_vely = vely(p.I);
                my_velz = velz(p.I);
              } else {
                if (tov_lapse)
                  alp_cell(p.I) = exp(phi_point[star]);
                if (tov_shift) {
                  betax_cell(p.I) = 0.0;
                  betay_cell(p.I) = 0.0;
                  betaz_cell(p.I) = 0.0;
                }
                my_velx = TOV_Velocity_x[star];
                my_vely = TOV_Velocity_y[star];
                my_velz = TOV_Velocity_z[star];
              }
            } else {
              if (tov_lapse)
                alp_cell(p.I) = exp(phi_point[star]);
              if (tov_shift) {
                betax_cell(p.I) = 0.0;
                betay_cell(p.I) = 0.0;
                betaz_cell(p.I) = 0.0;
              }
              my_velx = TOV_Velocity_x[star];
              my_vely = TOV_Velocity_y[star];
              my_velz = TOV_Velocity_z[star];
            }
          } else /* do not use old initial data */
          {
            gxx_cell(p.I) = max_g_diag;
            gyy_cell(p.I) = max_g_diag;
            gzz_cell(p.I) = max_g_diag;
            gxy_cell(p.I) = gxz_cell(p.I) = gyz_cell(p.I) = 0.0;
            if (tov_lapse)
              alp_cell(p.I) = exp(phi_point[star]);
            if (tov_shift) {
              betax_cell(p.I) = 0.0;
              betay_cell(p.I) = 0.0;
              betaz_cell(p.I) = 0.0;
            }
            my_velx = TOV_Velocity_x[star];
            my_vely = TOV_Velocity_y[star];
            my_velz = TOV_Velocity_z[star];
          }
          /* set to defined velocity. default is 0.0 because other velocities
           * violate Einsteins equations */
          velx(p.I) = my_velx;
          vely(p.I) = my_vely;
          velz(p.I) = my_velz;
          rho(p.I) = max_rho;
          eps(p.I) = eps_point[star];
          press(p.I) = press_point[star];

        } else if (combine_method == tov_combine_method_t::average) {
          /* here we 'average' all values in a more intelligent way */
          if (TOV_Use_Old_Matter_Initial_Data)
            max_rho = rho(p.I);
          else {
            max_rho = 0.0;
            rho(p.I) = 0.0;
          }
          star = -1;
          for (int star_i = 0; star_i < TOV_Num_TOVs; star_i++) {
            gxx_cell(p.I) +=
                (r_point[star_i] * r_point[star_i] /
                 (r_to_star[star_i] * r_to_star[star_i] + 1.0e-30)) /
                my_psi4;
            if (tov_lapse)
              alp_cell(p.I) *= exp(phi_point[star_i]);
            if (tov_shift) {
              betax_cell(p.I) = 0.0;
              betay_cell(p.I) = 0.0;
              betaz_cell(p.I) = 0.0;
            }
            rho(p.I) += rho_point[star_i];
            eps(p.I) += eps_point[star_i];
            press(p.I) += press_point[star_i];
            /* we still have to know if we are inside one star - and which */
            if (rho_point[star_i] > max_rho) {
              max_rho = rho_point[star_i];
              star = star_i;
            }
          }
          gxx_cell(p.I) -=
              ((TOV_Num_TOVs + TOV_Use_Old_Initial_Data - 1) / my_psi4);
          gyy_cell(p.I) = gxx_cell(p.I);
          gzz_cell(p.I) = gxx_cell(p.I);
          /* set to defined velocity. default is 0.0 because other velocities
           * violate the constraints */
          if (star > -1) {
            velx(p.I) = TOV_Velocity_x[star];
            vely(p.I) = TOV_Velocity_y[star];
            velz(p.I) = TOV_Velocity_z[star];
          }
        }
      });

  /* free local arrays */
  free(r_to_star);
  free(press_point);
  free(rho_point);
  free(eps_point);
  free(mu_point);
  free(phi_point);
  free(r_point);
} // TOVX_C_Exact

extern "C" void TOVX_C_Exact_Interpolation_C2V(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TOVX_C_Exact_Interpolation_C2V;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<0, 0, 0>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               gxx(p.I) = calc_avg_c2v(gxx_cell, p);
                               gxy(p.I) = calc_avg_c2v(gxy_cell, p);
                               gxz(p.I) = calc_avg_c2v(gxz_cell, p);
                               gyy(p.I) = calc_avg_c2v(gyy_cell, p);
                               gyz(p.I) = calc_avg_c2v(gyz_cell, p);
                               gzz(p.I) = calc_avg_c2v(gzz_cell, p);

                               kxx(p.I) = calc_avg_c2v(kxx_cell, p);
                               kxy(p.I) = calc_avg_c2v(kxy_cell, p);
                               kxz(p.I) = calc_avg_c2v(kxz_cell, p);
                               kyy(p.I) = calc_avg_c2v(kyy_cell, p);
                               kyz(p.I) = calc_avg_c2v(kyz_cell, p);
                               kzz(p.I) = calc_avg_c2v(kzz_cell, p);

                               alp(p.I) = calc_avg_c2v(alp_cell, p);
                               betax(p.I) = calc_avg_c2v(betax_cell, p);
                               betay(p.I) = calc_avg_c2v(betay_cell, p);
                               betaz(p.I) = calc_avg_c2v(betaz_cell, p);

                               dtalp(p.I) = calc_avg_c2v(dtalp_cell, p);
                               dtbetax(p.I) = calc_avg_c2v(dtbetax_cell, p);
                               dtbetay(p.I) = calc_avg_c2v(dtbetay_cell, p);
                               dtbetaz(p.I) = calc_avg_c2v(dtbetaz_cell, p);
                             });

  CCTK_INFO("Done interpolation for TOV hydro initial data.");

} // TOVX_C_Exact_Interpolation_C2V

/*@@
   @routine    TOVX_Exact_ADM
   @date       Thu Oct 24 14:30:00 2002
   @author     Frank Loeffler, converted fortran routine by Ian Hawke
   @desc
       Schedule routine for interpolation of 1D to 3D grid for only ADM
variables
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
@@*/
extern "C" void TOVX_C_Exact_ADM(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TOVX_C_Exact_ADM;
  DECLARE_CCTK_PARAMETERS;

  assert(TOV_Surface != 0);
  assert(TOV_R_Surface != 0);
  assert(TOV_r_1d != 0);
  assert(TOV_rbar_1d != 0);
  assert(TOV_press_1d != 0);
  assert(TOV_phi_1d != 0);
  assert(TOV_m_1d != 0);

  /* allocate local arrays */
  CCTK_REAL *r_to_star = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));
  CCTK_REAL *press_point = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));
  CCTK_REAL *rho_point = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));
  CCTK_REAL *eps_point = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));
  CCTK_REAL *mu_point = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));
  CCTK_REAL *phi_point = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));
  CCTK_REAL *r_point = (CCTK_REAL *)calloc(TOV_Num_TOVs, sizeof(CCTK_REAL));

  /* clear initial data */
  if (TOV_Clear_Initial_Data > 0 && !(TOV_Use_Old_Initial_Data)) {
    grid.loop_int<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          alp_cell(p.I) = 1.0;
          betax_cell(p.I) = 0.0;
          betay_cell(p.I) = 0.0;
          betaz_cell(p.I) = 0.0;

          dtalp_cell(p.I) = 0.0;
          dtbetax_cell(p.I) = 0.0;
          dtbetay_cell(p.I) = 0.0;
          dtbetaz_cell(p.I) = 0.0;

          gxx_cell(p.I) = 0.0;
          gyy_cell(p.I) = 0.0;
          gzz_cell(p.I) = 0.0;
          gxy_cell(p.I) = gxz_cell(p.I) = gyz_cell(p.I) = 0.0;

          kxx_cell(p.I) = 0.0;
          kyy_cell(p.I) = 0.0;
          kzz_cell(p.I) = 0.0;
          kxy_cell(p.I) = kxz_cell(p.I) = kyz_cell(p.I) = 0.0;
        });
  }

  CCTK_INT tov_lapse = CCTK_EQUALS(initial_lapse, "tov");
  CCTK_INT tov_shift = CCTK_EQUALS(initial_shift, "tov");

  tov_combine_method_t combine_method;
  if (CCTK_EQUALS(TOV_Combine_Method, "maximum")) {
    combine_method = tov_combine_method_t::maximum;
  } else if (CCTK_EQUALS(TOV_Combine_Method, "average")) {
    combine_method = tov_combine_method_t::average;
  } else {
    CCTK_ERROR("Unknown value for parameter \"TOV_Combine_Method\"");
  }

  auto TOV_press_1d_local = TOV_press_1d;
  auto TOV_phi_1d_local = TOV_phi_1d;
  auto TOV_rbar_1d_local = TOV_rbar_1d;
  auto TOV_r_1d_local = TOV_r_1d;
  auto TOV_Surface_local = TOV_Surface;
  auto TOV_R_Surface_local = TOV_R_Surface;

  grid.loop_int<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_INT star;
        CCTK_REAL g_diag, max_g_diag, max_rho;
        CCTK_REAL local_tiny = 1.0e-14;
        /* remember the old conformal factor to the power of 4 */
        CCTK_REAL my_psi4 = 1.0;

        for (star = 0; star < TOV_Num_TOVs; star++) {
          r_to_star[star] =
              sqrt((p.x - TOV_Position_x[star]) * (p.x - TOV_Position_x[star]) +
                   (p.y - TOV_Position_y[star]) * (p.y - TOV_Position_y[star]) +
                   (p.z - TOV_Position_z[star]) * (p.z - TOV_Position_z[star]));
          int star_i = star * TOV_Num_Radial;

          /* do the actual interpolation */
          TOVX_C_interp_tov_isotropic(
              star, TOV_Num_Radial, TOV_Fast_Interpolation,
              &(TOV_press_1d_local[star_i]), &(TOV_phi_1d_local[star_i]),
              &(TOV_rbar_1d_local[star_i]), &(TOV_r_1d_local[star_i]),
              &(r_to_star[star]), TOV_Surface_local[star], &(press_point[star]),
              &(phi_point[star]), &(r_point[star]));

          /* is some perturbation wanted? */
          if (Perturb[star] == 0)
            rho_point[star] = pow(press_point[star] / TOV_K, 1.0 / TOV_Gamma);
          else {
            CCTK_REAL r = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
            CCTK_REAL PI = 4.0 * atan(1.0);
            rho_point[star] =
                pow(press_point[star] / TOV_K, 1.0 / TOV_Gamma) *
                (1.0 + Pert_Amplitude[star] *
                           cos(PI / 2.0 * r / TOV_R_Surface_local[star]));
          }
          /* set eps_point */
          if (rho_point[star] > local_tiny)
            eps_point[star] =
                press_point[star] / (TOV_Gamma - 1.0) / rho_point[star];
          else
            eps_point[star] = 0.0;
          /* set mu_point */
          mu_point[star] = rho_point[star] * (1.0 + eps_point[star]);
        }

        /* find out from which star we want to have the data */
        if (combine_method == tov_combine_method_t::maximum) {
          /* to do this, we use here simply the max of the gxx-value */
          star = 0;
          max_g_diag = 0.0;
          max_rho = rho_point[0];
          for (int star_i = 0; star_i < TOV_Num_TOVs; star_i++) {
            g_diag = (r_point[star_i] / (r_to_star[star_i] + 1.0e-30)) *
                     (r_point[star_i] / (r_to_star[star_i] + 1.0e-30));
            if ((g_diag - max_g_diag) > local_tiny) {
              max_g_diag = g_diag;
              star = star_i;
            }
            if ((rho_point[star_i] - max_rho) > local_tiny) {
              max_rho = rho_point[star_i];
              star = star_i;
            }
          }

          /* handle initial data */
          if (TOV_Use_Old_Initial_Data) {
            /* check metric */
            if ((my_psi4 * gxx_cell(p.I) < max_g_diag) &&
                (my_psi4 * gyy_cell(p.I) < max_g_diag) &&
                (my_psi4 * gzz_cell(p.I) < max_g_diag)) {
              gxx_cell(p.I) = max_g_diag / my_psi4;
              gyy_cell(p.I) = max_g_diag / my_psi4;
              gzz_cell(p.I) = max_g_diag / my_psi4;
              gxy_cell(p.I) = gxz_cell(p.I) = gyz_cell(p.I) = 0.0;
            }
            if (tov_lapse)
              alp_cell(p.I) = exp(phi_point[star]);
            if (tov_shift) {
              betax_cell(p.I) = 0.0;
              betay_cell(p.I) = 0.0;
              betaz_cell(p.I) = 0.0;
            }
          } else /* do not use old initial data */
          {
            /* no psi, since it is 1.0 here */
            /* but maybe we want to have it != 1.0 */
            gxx_cell(p.I) = max_g_diag;
            gyy_cell(p.I) = max_g_diag;
            gzz_cell(p.I) = max_g_diag;
            gxy_cell(p.I) = gxz_cell(p.I) = gyz_cell(p.I) = 0.0;
            if (tov_lapse)
              alp_cell(p.I) = exp(phi_point[star]);
            if (tov_shift) {
              betax_cell(p.I) = 0.0;
              betay_cell(p.I) = 0.0;
              betaz_cell(p.I) = 0.0;
            }
          }

        } else if (combine_method == tov_combine_method_t::average) {
          /* here we 'average' all values in a more intelligent way */
          star = -1;
          for (int star_i = 0; star_i < TOV_Num_TOVs; star_i++) {
            /* Conformal_Flat_Three_Metric unsupported yet */
            gxx_cell(p.I) +=
                (r_point[star_i] * r_point[star_i] /
                 (r_to_star[star_i] * r_to_star[star_i] + 1.0e-30)) /
                my_psi4;
            if (tov_lapse)
              alp_cell(p.I) *= exp(phi_point[star_i]);
            if (tov_shift) {
              betax_cell(p.I) = 0.0;
              betay_cell(p.I) = 0.0;
              betaz_cell(p.I) = 0.0;
            }
          }
          gxx_cell(p.I) -=
              ((TOV_Num_TOVs + TOV_Use_Old_Initial_Data - 1) / my_psi4);
          gyy_cell(p.I) = gxx_cell(p.I);
          gzz_cell(p.I) = gxx_cell(p.I);
        }
      });

  /* free local arrays */
  free(r_to_star);
  free(press_point);
  free(rho_point);
  free(eps_point);
  free(mu_point);
  free(phi_point);
  free(r_point);
} // TOV_Exact_ADM

extern "C" void TOVX_C_Exact_ADM_Interpolation_C2V(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TOVX_C_Exact_ADM_Interpolation_C2V;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<0, 0, 0>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               gxx(p.I) = calc_avg_c2v(gxx_cell, p);
                               gxy(p.I) = calc_avg_c2v(gxy_cell, p);
                               gxz(p.I) = calc_avg_c2v(gxz_cell, p);
                               gyy(p.I) = calc_avg_c2v(gyy_cell, p);
                               gyz(p.I) = calc_avg_c2v(gyz_cell, p);
                               gzz(p.I) = calc_avg_c2v(gzz_cell, p);

                               kxx(p.I) = calc_avg_c2v(kxx_cell, p);
                               kxy(p.I) = calc_avg_c2v(kxy_cell, p);
                               kxz(p.I) = calc_avg_c2v(kxz_cell, p);
                               kyy(p.I) = calc_avg_c2v(kyy_cell, p);
                               kyz(p.I) = calc_avg_c2v(kyz_cell, p);
                               kzz(p.I) = calc_avg_c2v(kzz_cell, p);

                               alp(p.I) = calc_avg_c2v(alp_cell, p);
                               betax(p.I) = calc_avg_c2v(betax_cell, p);
                               betay(p.I) = calc_avg_c2v(betay_cell, p);
                               betaz(p.I) = calc_avg_c2v(betaz_cell, p);

                               dtalp(p.I) = calc_avg_c2v(dtalp_cell, p);
                               dtbetax(p.I) = calc_avg_c2v(dtbetax_cell, p);
                               dtbetay(p.I) = calc_avg_c2v(dtbetay_cell, p);
                               dtbetaz(p.I) = calc_avg_c2v(dtbetaz_cell, p);
                             });

  CCTK_INFO("Done interpolation for TOV ADM variables.");

} // TOV_Exact_ADM_Interpolation_C2V

/*@@
   @routine    TOVX_Integrate_RHS
   @date       Thu Oct 24 14:30:00 2002
   @author     Frank Loeffler, converted fortran routine by Ian Hawke
   @desc
      Integrates the ODEs using RK4.
      We rescale at the end to match to a Schwarzschild exterior.
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
@@*/
extern "C" void TOVX_C_Integrate_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TOVX_C_Integrate_RHS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT star, star_i, i, TOV_Surface_Index;
  CCTK_REAL old_data[NUMVARS], source_data[NUMVARS], in_data[NUMVARS],
      new_data[NUMVARS], k1[NUMVARS], k2[NUMVARS], k3[NUMVARS], k4[NUMVARS];
  CCTK_REAL Surface_Mass, factor, local_rho;
  CCTK_REAL LOCAL_TINY = 1.0e-20;

  assert(TOV_Surface != 0);
  assert(TOV_R_Surface != 0);
  assert(TOV_RProp_Surface != 0);

  assert(TOV_r_1d != 0);
  assert(TOV_rbar_1d != 0);
  assert(TOV_press_1d != 0);
  assert(TOV_phi_1d != 0);
  assert(TOV_m_1d != 0);
  assert(TOV_mbary_1d != 0);
  assert(TOV_rprop_1d != 0);

  /* do it for all stars */
  for (star = 0; star < TOV_Num_TOVs; star++) {
    /* remember array index */
    star_i = star * TOV_Num_Radial;
    const CCTK_REAL rho_central = TOV_Rho_Central[star];

    /* Set conformal state like set in parameter file if we do not use
     * the old initial data. In this case we have to use what we get */

    /* clear arrays first */
    TOVX_C_fill(&(TOV_press_1d[star_i]), TOV_Num_Radial, 0.0);
    TOVX_C_fill(&(TOV_m_1d[star_i]), TOV_Num_Radial, 0.0);
    TOVX_C_fill(&(TOV_phi_1d[star_i]), TOV_Num_Radial, 0.0);
    TOVX_C_fill(&(TOV_rbar_1d[star_i]), TOV_Num_Radial, 0.0);
    TOVX_C_fill(&(TOV_r_1d[star_i]), TOV_Num_Radial, 0.0);
    TOVX_C_fill(&(TOV_mbary_1d[star_i]), TOV_Num_Radial, 0.0);
    TOVX_C_fill(&(TOV_rprop_1d[star_i]), TOV_Num_Radial, 0.0);

    /* set start values */
    TOV_press_1d[star_i] = TOV_K * pow(rho_central, TOV_Gamma);
    /* TOV_r_1d    [star_i] = LOCAL_TINY;
    TOV_rbar_1d [star_i] = LOCAL_TINY;*/

    /* build TOV_r_1d[] */
    for (i = star_i + 1; i < star_i + TOV_Num_Radial; i++)
      TOV_r_1d[i] = TOV_r_1d[i - 1] + TOV_dr[star];

    TOV_Surface[star] = -1.0;
    TOV_Surface_Index = -1.0;

#define RKLOOP for (int rk = 0; rk < NUMVARS; rk++)
    /* loop over all radii */
    for (i = star_i;
         (i < star_i + TOV_Num_Radial - 1) && (TOV_Surface[star] < 0.0); i++) {
      /* set up RK arrays */
      old_data[0] = TOV_press_1d[i];
      old_data[1] = TOV_m_1d[i];
      old_data[2] = TOV_phi_1d[i];
      if (fabs(TOV_rbar_1d[i] - TOV_r_1d[i]) < LOCAL_TINY)
        old_data[3] = 0.0;
      else
        old_data[3] = log(TOV_rbar_1d[i] / TOV_r_1d[i]);
      old_data[4] = TOV_mbary_1d[i];
      old_data[5] = TOV_rprop_1d[i];

      /* usual RK4 */
      RKLOOP in_data[rk] = old_data[rk];
      TOVX_C_fill(source_data, 6, 0.0);

      TOVX_C_Source_RHS(TOV_r_1d[i], TOV_K, TOV_Gamma, in_data, source_data);

      RKLOOP k1[rk] = TOV_dr[star] * source_data[rk];
      RKLOOP in_data[rk] = old_data[rk] + 0.5 * k1[rk];
      TOVX_C_Source_RHS(TOV_r_1d[i] + 0.5 * TOV_dr[star], TOV_K, TOV_Gamma,
                        in_data, source_data);

      RKLOOP k2[rk] = TOV_dr[star] * source_data[rk];
      RKLOOP in_data[rk] = old_data[rk] + 0.5 * k2[rk];
      TOVX_C_Source_RHS(TOV_r_1d[i] + 0.5 * TOV_dr[star], TOV_K, TOV_Gamma,
                        in_data, source_data);

      RKLOOP k3[rk] = TOV_dr[star] * source_data[rk];
      RKLOOP in_data[rk] = old_data[rk] + k3[rk];
      TOVX_C_Source_RHS(TOV_r_1d[i] + TOV_dr[star], TOV_K, TOV_Gamma, in_data,
                        source_data);
      RKLOOP k4[rk] = TOV_dr[star] * source_data[rk];
      RKLOOP new_data[rk] =
          old_data[rk] + (k1[rk] + k4[rk] + 2.0 * (k2[rk] + k3[rk])) / 6.0;

      TOV_press_1d[i + 1] = new_data[0];
      TOV_m_1d[i + 1] = new_data[1];
      TOV_phi_1d[i + 1] = new_data[2];
      TOV_rbar_1d[i + 1] = TOV_r_1d[i + 1] * exp(new_data[3]);
      TOV_mbary_1d[i + 1] = new_data[4];
      TOV_rprop_1d[i + 1] = new_data[5];

      /* otherwise the code crashes later */
      if (TOV_press_1d[i + 1] < 0.0)
        TOV_press_1d[i + 1] = 0.0;

      local_rho = pow(TOV_press_1d[i + 1] / TOV_K, 1.0 / TOV_Gamma);

      /* scan for the surface */
      if ((local_rho <= 0.0) || (TOV_press_1d[i + 1] <= 0.0)) {
        TOV_Surface[star] = TOV_r_1d[i];
        TOV_R_Surface[star] = TOV_rbar_1d[i];
        TOV_RProp_Surface[star] = TOV_rprop_1d[i];
        TOV_Surface_Index = i;
      }
    }
    if (TOV_Surface[star] < 0.0)
      CCTK_WARN(0, "Did not integrate out to surface of the star! "
                   "Increase TOV_dr or TOV_Num_Radial and rerun");

    Surface_Mass = TOV_m_1d[TOV_Surface_Index];
    factor =
        0.5 *
        (sqrt(TOV_Surface[star] * (TOV_Surface[star] - 2.00 * Surface_Mass)) +
         TOV_Surface[star] - Surface_Mass) /
        TOV_rbar_1d[TOV_Surface_Index];

    TOV_R_Surface[star] *= factor;
    for (i = star_i; i < star_i + TOV_Num_Radial; i++) {
      TOV_rbar_1d[i] *= factor;
      TOV_phi_1d[i] -= TOV_phi_1d[TOV_Surface_Index] -
                       0.5 * log(1.0 - 2.0 * Surface_Mass / TOV_Surface[star]);
      /* match to Schwarzschield */
      if (i > TOV_Surface_Index) {
        TOV_press_1d[i] = 0.0;
        TOV_rbar_1d[i] =
            0.5 * (sqrt(TOV_r_1d[i] * (TOV_r_1d[i] - 2.0 * Surface_Mass)) +
                   TOV_r_1d[i] - Surface_Mass);
        TOV_m_1d[i] = Surface_Mass;
        TOV_phi_1d[i] = 0.5 * log(1.0 - 2.0 * Surface_Mass / TOV_r_1d[i]);
        TOV_mbary_1d[i] = TOV_mbary_1d[TOV_Surface_Index];
      }
    }
  }
  CCTK_INFO("Integrated TOV equation");
  /* do some info */
  CCTK_VInfo(CCTK_THORNSTRING, "Information about the TOVs used:");
  CCTK_VInfo("", "TOV    radius    mass  bary_mass mass(g) cent.rho rho(cgi)   "
                 "     K   K(cgi)    Gamma");
  for (i = 0; i < TOV_Num_TOVs; i++)
    if (fabs(TOV_Gamma - 2.0) < LOCAL_TINY)
      CCTK_VInfo("", "  %d  %8g %8g %8g %8.3g %8g %8.3g %8g %8.3g %8g",
                 (int)i + 1, TOV_R_Surface[i],
                 TOV_m_1d[(i + 1) * TOV_Num_Radial - 1],
                 TOV_mbary_1d[(i + 1) * TOV_Num_Radial - 1],
                 TOV_m_1d[(i + 1) * TOV_Num_Radial - 1] * CONSTANT_Msolar_cgi,
                 TOV_Rho_Central[i],
                 TOV_Rho_Central[i] / pow(CONSTANT_G_cgi, 3.0) /
                     pow(CONSTANT_Msolar_cgi, 2.0) * pow(CONSTANT_c_cgi, 6.0),
                 TOV_K,
                 TOV_K * pow(CONSTANT_G_cgi, 3.0) *
                     pow(CONSTANT_Msolar_cgi, 2.0) / pow(CONSTANT_c_cgi, 4.0),
                 TOV_Gamma);
    else
      CCTK_VInfo("", "  %d  %8g %8g %8.3g %8g %8.3g %8g %8g", (int)i + 1,
                 TOV_R_Surface[i], TOV_m_1d[(i + 1) * TOV_Num_Radial - 1],
                 TOV_m_1d[(i + 1) * TOV_Num_Radial - 1] * CONSTANT_Msolar_cgi,
                 TOV_Rho_Central[i],
                 TOV_Rho_Central[i] / pow(CONSTANT_G_cgi, 3.0) /
                     pow(CONSTANT_Msolar_cgi, 2.0) * pow(CONSTANT_c_cgi, 6.0),
                 TOV_K, TOV_Gamma);

  CCTK_INFO("Done Integrating for TOV 1D data");
} // TOVX_C_Integrate_RHS

extern "C" void TOVX_C_AllocateMemory(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TOVX_C_AllocateMemory;
  DECLARE_CCTK_PARAMETERS;

  assert(TOV_Surface == 0);
  assert(TOV_R_Surface == 0);
  assert(TOV_RProp_Surface == 0);

  assert(TOV_r_1d == 0);
  assert(TOV_rbar_1d == 0);
  assert(TOV_press_1d == 0);
  assert(TOV_phi_1d == 0);
  assert(TOV_m_1d == 0);
  assert(TOV_mbary_1d == 0);
  assert(TOV_rprop_1d == 0);

  TOV_Surface =
      (CCTK_REAL *)amrex::The_Arena()->alloc(TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_R_Surface =
      (CCTK_REAL *)amrex::The_Arena()->alloc(TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_RProp_Surface =
      (CCTK_REAL *)amrex::The_Arena()->alloc(TOV_Num_TOVs * sizeof(CCTK_REAL));

  TOV_r_1d = (CCTK_REAL *)amrex::The_Arena()->alloc(
      TOV_Num_Radial * TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_rbar_1d = (CCTK_REAL *)amrex::The_Arena()->alloc(
      TOV_Num_Radial * TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_press_1d = (CCTK_REAL *)amrex::The_Arena()->alloc(
      TOV_Num_Radial * TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_phi_1d = (CCTK_REAL *)amrex::The_Arena()->alloc(
      TOV_Num_Radial * TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_m_1d = (CCTK_REAL *)amrex::The_Arena()->alloc(
      TOV_Num_Radial * TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_mbary_1d = (CCTK_REAL *)amrex::The_Arena()->alloc(
      TOV_Num_Radial * TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_rprop_1d = (CCTK_REAL *)amrex::The_Arena()->alloc(
      TOV_Num_Radial * TOV_Num_TOVs * sizeof(CCTK_REAL));
} // TOVX_C_AllocateMemory

extern "C" void TOVX_C_FreeMemory(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TOVX_C_AllocateMemory;
  DECLARE_CCTK_PARAMETERS;

  assert(TOV_Surface != 0);
  assert(TOV_R_Surface != 0);
  assert(TOV_RProp_Surface != 0);

  assert(TOV_r_1d != 0);
  assert(TOV_rbar_1d != 0);
  assert(TOV_press_1d != 0);
  assert(TOV_phi_1d != 0);
  assert(TOV_m_1d != 0);

  amrex::The_Arena()->free(TOV_Surface);
  amrex::The_Arena()->free(TOV_R_Surface);
  amrex::The_Arena()->free(TOV_RProp_Surface);

  amrex::The_Arena()->free(TOV_r_1d);
  amrex::The_Arena()->free(TOV_rbar_1d);
  amrex::The_Arena()->free(TOV_press_1d);
  amrex::The_Arena()->free(TOV_phi_1d);
  amrex::The_Arena()->free(TOV_m_1d);
  amrex::The_Arena()->free(TOV_mbary_1d);
  amrex::The_Arena()->free(TOV_rprop_1d);

  TOV_Surface = 0;
  TOV_R_Surface = 0;
  TOV_RProp_Surface = 0;

  TOV_r_1d = 0;
  TOV_rbar_1d = 0;
  TOV_press_1d = 0;
  TOV_phi_1d = 0;
  TOV_m_1d = 0;
  TOV_mbary_1d = 0;
  TOV_rprop_1d = 0;
} // TOVX_FreeMemory

#include "external.hxx"
} // namespace TOVSolverX
