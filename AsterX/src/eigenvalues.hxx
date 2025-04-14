#ifndef ASTERX_EIGENVALUES_HXX
#define ASTERX_EIGENVALUES_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

namespace AsterX {
using namespace std;
using namespace Arith;

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec<CCTK_REAL, 3>
compute_a(CCTK_REAL alp_avg, CCTK_REAL beta_avg, CCTK_REAL u_avg, CCTK_REAL vel,
          CCTK_REAL rho, CCTK_REAL cs2, CCTK_REAL w_lor, CCTK_REAL h,
          CCTK_REAL bsq) {

  const CCTK_REAL vA2 = bsq / (rho * h + bsq);
  const CCTK_REAL v02 = vA2 + cs2 * (1 - vA2);
  const CCTK_REAL u02 = pow2(w_lor / alp_avg);
  const CCTK_REAL one_over_alp2 = 1.0 / pow2(alp_avg);

  vec<CCTK_REAL, 3> a{
      (1.0 - v02) * u02 + v02 * one_over_alp2,
      2.0 * (beta_avg * one_over_alp2 * v02 - u02 * vel * (1.0 - v02)),
      u02 * pow2(vel) * (1.0 - v02) -
          v02 * (u_avg - pow2(beta_avg) * one_over_alp2)};

  return a;
}

inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST vec<vec<CCTK_REAL, 4>, 2>
    eigenvalues(CCTK_REAL alp_avg, CCTK_REAL beta_avg, CCTK_REAL u_avg,
                vec<CCTK_REAL, 2> vel, vec<CCTK_REAL, 2> rho,
                vec<CCTK_REAL, 2> cs2, vec<CCTK_REAL, 2> w_lor,
                vec<CCTK_REAL, 2> h, vec<CCTK_REAL, 2> bsq) {
  // computing characteristics for the minus side
  // See Eq. (28) of Giacomazzo & Rezzolla (2007) with b^i=0

  vec<CCTK_REAL, 3> a_m = compute_a(alp_avg, beta_avg, u_avg, vel(0), rho(0),
                                    cs2(0), w_lor(0), h(0), bsq(0));

  CCTK_REAL det_m = pow2(a_m(1)) - 4 * a_m(2) * a_m(0);
  CCTK_REAL sqrt_det_m = sqrt(0.5 * (det_m + fabs(det_m)));

  vec<CCTK_REAL, 4> lambda_m{((-a_m(1) + sqrt_det_m) / (2 * a_m(0))),
                             ((-a_m(1) + sqrt_det_m) / (2 * a_m(0))),
                             ((-a_m(1) - sqrt_det_m) / (2 * a_m(0))),
                             ((-a_m(1) - sqrt_det_m) / (2 * a_m(0)))};

  // computing characteristics for the plus side

  vec<CCTK_REAL, 3> a_p = compute_a(alp_avg, beta_avg, u_avg, vel(1), rho(1),
                                    cs2(1), w_lor(1), h(1), bsq(1));

  CCTK_REAL det_p = pow2(a_p(1)) - 4 * a_p(2) * a_p(0);
  CCTK_REAL sqrt_det_p = sqrt(0.5 * (det_p + fabs(det_p)));

  vec<CCTK_REAL, 4> lambda_p{((-a_p(1) + sqrt_det_p) / (2 * a_p(0))),
                             ((-a_p(1) + sqrt_det_p) / (2 * a_p(0))),
                             ((-a_p(1) - sqrt_det_p) / (2 * a_p(0))),
                             ((-a_p(1) - sqrt_det_p) / (2 * a_p(0)))};

  // 2D array containing characteristics for left (minus) and right
  // (plus) sides
  vec<vec<CCTK_REAL, 4>, 2> lambda{lambda_m, lambda_p};
  return lambda;
};

} // namespace AsterX

#endif // ASTERX_EIGENVALUES_HXX
