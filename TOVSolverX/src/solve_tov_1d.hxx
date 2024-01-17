/* file    solve_tov_1d.hxx
 * author  Frank Loeffler, converted from fortran thorn by Ian Hawk
 *                         modified to be compatible with CarpetX by Johnny
 *                         original thorn from
 * Cactus/arrangements/EinsteinInitialData/TOVSolverX date    2022/07/24 desc
 * Helper functions for solving 1D TOV equations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "constants.h"
#include "utils.hxx"
#define NUMVARS 6

/* centered differencing with one-sided differencing at the boundary */
#define DIFF_X(a)                                                              \
  (((i == 0) ? (a[CCTK_GFINDEX3D(cctkGH, i + 1, j, k)] -                       \
                a[CCTK_GFINDEX3D(cctkGH, i, j, k)])                            \
    : (i == (cctk_lsh[0] - 1))                                                 \
        ? (a[CCTK_GFINDEX3D(cctkGH, i, j, k)] -                                \
           a[CCTK_GFINDEX3D(cctkGH, i - 1, j, k)])                             \
        : 0.5 * (a[CCTK_GFINDEX3D(cctkGH, i + 1, j, k)] -                      \
                 a[CCTK_GFINDEX3D(cctkGH, i - 1, j, k)])) /                    \
   CCTK_DELTA_SPACE(0))
#define DIFF_Y(a)                                                              \
  (((j == 0) ? (a[CCTK_GFINDEX3D(cctkGH, i, j + 1, k)] -                       \
                a[CCTK_GFINDEX3D(cctkGH, i, j, k)])                            \
    : (j == (cctk_lsh[1] - 1))                                                 \
        ? (a[CCTK_GFINDEX3D(cctkGH, i, j, k)] -                                \
           a[CCTK_GFINDEX3D(cctkGH, i, j - 1, k)])                             \
        : 0.5 * (a[CCTK_GFINDEX3D(cctkGH, i, j + 1, k)] -                      \
                 a[CCTK_GFINDEX3D(cctkGH, i, j - 1, k)])) /                    \
   CCTK_DELTA_SPACE(1))
#define DIFF_Z(a)                                                              \
  (((k == 0) ? (a[CCTK_GFINDEX3D(cctkGH, i, j, k + 1)] -                       \
                a[CCTK_GFINDEX3D(cctkGH, i, j, k)])                            \
    : (k == (cctk_lsh[2] - 1))                                                 \
        ? (a[CCTK_GFINDEX3D(cctkGH, i, j, k)] -                                \
           a[CCTK_GFINDEX3D(cctkGH, i, j, k - 1)])                             \
        : 0.5 * (a[CCTK_GFINDEX3D(cctkGH, i, j, k + 1)] -                      \
                 a[CCTK_GFINDEX3D(cctkGH, i, j, k - 1)])) /                    \
   CCTK_DELTA_SPACE(2))

namespace TOVSolverX {

/*@@
   @routine    TOVX_Source_RHS
   @date       Thu Oct 24 14:30:00 2002
   @author     Frank Loeffler - converted fortran routine by Ian Hawke
   @desc
      The source terms for the ODEs. These are equations (2), (3), (4)
      and (18) from the Baumgarte notes.
      That is the vector in order is (P, m, phi, rbar).
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
@@*/
CCTK_HOST void TOVX_C_Source_RHS(CCTK_REAL r, CCTK_REAL K, CCTK_REAL Gamma,
                                CCTK_REAL old_data[NUMVARS],
                                CCTK_REAL source_data[NUMVARS]) {
  CCTK_REAL LOCAL_TINY, PI;
  CCTK_REAL press, rho, eps, mu, m;
  CCTK_REAL r_minus_two_m;

  LOCAL_TINY = 1.0e-35;
  PI = 4.0 * atan(1.0);

  press = old_data[0];
  if (press < LOCAL_TINY)
    press = LOCAL_TINY;
  m = old_data[1];

  rho = pow(press / K, 1.0 / Gamma);
  eps = press / (Gamma - 1.0) / rho;
  mu = rho * (1.0 + eps);

  r_minus_two_m = r - 2.0 * m;

  if ((r <= 0.0) && (m <= 0.0)) {
    source_data[1] = 0.0;
    source_data[2] = 0.0;
    source_data[3] = 0.0;
    source_data[4] = 0.0;
    source_data[5] = 0.0;
  } else {
    source_data[2] = (m + 4 * PI * r * r * r * press) / r_minus_two_m / r;
    /* source_data[0] = -(press + mu) * source_data[2]; */
    source_data[0] =
        -(press + mu) * (m + 4 * PI * r * r * r * press) / r_minus_two_m / r;
    source_data[1] = 4 * PI * r * r * mu;
    source_data[3] = (sqrt(r) - sqrt(r_minus_two_m)) / r / sqrt(r_minus_two_m);
    source_data[5] = 1.0 / sqrt(1.0 - 2.0 * m / r);
    source_data[4] = source_data[5] * 4 * PI * rho * r * r;
  }
}

/*----------------------------------------------------------------------------*/

/* utility routine
 * recursive search-routine for arrays
 * here used to look for the last index in an ordered array with its
 * value < goal
 */
CCTK_HOST CCTK_INT TOVX_C_find_index(CCTK_INT array_size, CCTK_REAL *array,
                                    CCTK_REAL goal, CCTK_INT lower_index,
                                    CCTK_INT upper_index) {
  CCTK_INT middle_index;

  if (lower_index >= (upper_index - 1))
    return lower_index;

  middle_index = (lower_index + upper_index) / 2;

  if (array[middle_index] < goal)
    return TOVX_C_find_index(array_size, array, goal, middle_index, upper_index);
  else
    return TOVX_C_find_index(array_size, array, goal, lower_index, middle_index);
}

/* utility rountine
 * interpolates from (thorn-internal) 1D-data to Cactus 3D-grid */
/* input is all but *press_point *phi_point and *r_point */
CCTK_HOST void TOVX_C_interp_tov_isotropic(
    CCTK_INT star, CCTK_INT TOV_Num_Radial, CCTK_INT TOV_Fast_Interpolation,
    CCTK_REAL *TOV_press_1d_local, CCTK_REAL *TOV_phi_1d_local,
    CCTK_REAL *TOV_rbar_1d_local, CCTK_REAL *TOV_r_1d_local, CCTK_REAL *r,
    CCTK_REAL surface, CCTK_REAL *press_point, CCTK_REAL *phi_point,
    CCTK_REAL *r_point) {
  CCTK_INT left_index;
  CCTK_REAL h, M;

  if (*r < 0.0)
    CCTK_WARN(0, "Negative radius found");
  if (*r < TOV_rbar_1d_local[1])
    *r = TOV_rbar_1d_local[1];
  if (*r > TOV_rbar_1d_local[TOV_Num_Radial - 2]) {
    {
      *press_point = 0.0;
      M = 0.5 * TOV_r_1d_local[TOV_Num_Radial - 1] *
          (1.0 - exp(2.0 * TOV_phi_1d_local[TOV_Num_Radial - 1]));
      *r_point = (2 * *r + M) * (2 * *r + M) * 0.25 / *r;
      *phi_point = 0.5 * log(1 - 2 * M / *r_point);
      return;
    }
  }

  if (TOV_Fast_Interpolation)
    left_index = TOVX_C_find_index(TOV_Num_Radial - 1, TOV_rbar_1d_local, *r, 0,
                                  TOV_Num_Radial - 1);
  else {
    left_index = 0;
    while ((left_index < TOV_Num_Radial - 2) &&
           (TOV_rbar_1d_local[left_index + 1] < *r))
      left_index++;
  }

  h = (*r - TOV_rbar_1d_local[left_index]) /
      (TOV_rbar_1d_local[left_index + 1] - TOV_rbar_1d_local[left_index]);
  *r_point = (1.0 - h) * TOV_r_1d_local[left_index] +
             h * TOV_r_1d_local[left_index + 1];
  *phi_point = (1.0 - h) * TOV_phi_1d_local[left_index] +
               h * TOV_phi_1d_local[left_index + 1];
  if (*r_point < surface)
    *press_point = (1.0 - h) * TOV_press_1d_local[left_index] +
                   h * TOV_press_1d_local[left_index + 1];
  else
    *press_point = 0.0;
}

} // namespace TOVSolverX
