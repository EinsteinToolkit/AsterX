#ifndef ASTERX_PRIM2CON_HXX
#define ASTERX_PRIM2CON_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>
#include "utils.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;

struct prim {
  CCTK_REAL rho;
  vec<CCTK_REAL, 3> vel;
  CCTK_REAL eps, press;
  vec<CCTK_REAL, 3> Bvec;
};

struct cons {
  CCTK_REAL dens;
  vec<CCTK_REAL, 3> mom;
  CCTK_REAL tau;
  vec<CCTK_REAL, 3> dBvec;
};

CCTK_DEVICE CCTK_HOST void prim2con(const smat<CCTK_REAL, 3> &g,
                                    const prim &pv, cons &cv) {

  // determinant of spatial metric
  const CCTK_REAL sqrt_detg = sqrt(calc_det(g));

  // TODO: compute specific internal energy based on user-specified EOS
  // currently, computing eps for classical ideal gas

  /* Computing v_j */
  const vec<CCTK_REAL, 3> &v_up = pv.vel;
  const vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);

  const CCTK_REAL w_lorentz = calc_wlorentz(v_low, v_up);

  /* Computing B_j */
  const vec<CCTK_REAL, 3> &B_up = pv.Bvec;
  const vec<CCTK_REAL, 3> B_low = calc_contraction(g, B_up);

  /* Computing b^t : this is b^0 * alp */
  const CCTK_REAL bst = w_lorentz * calc_contraction(B_up, v_low);

  /* Computing b_j */
  const vec<CCTK_REAL, 3> b_low = B_low / w_lorentz + bst * v_low;

  /* Computing b^mu b_mu */
  const CCTK_REAL bs2 =
      (calc_contraction(B_up, B_low) + bst * bst) / (w_lorentz * w_lorentz);

  // computing conserved from primitives
  cv.dens = sqrt_detg * pv.rho * w_lorentz;

  cv.mom = sqrt_detg * (w_lorentz * w_lorentz *
                            (pv.rho * (1.0 + pv.eps) + pv.press + bs2) * v_low -
                        bst * b_low);

  cv.tau = sqrt_detg * (w_lorentz * w_lorentz *
                            (pv.rho * (1.0 + pv.eps) + pv.press + bs2) -
                        (pv.press + 0.5 * bs2) - bst * bst) -
           cv.dens;

  cv.dBvec = sqrt_detg * pv.Bvec;
}

} // namespace AsterX

#endif // ASTERX_PRIM2CON_HXX
