#ifndef RECONX_EPPM_HXX
#define RECONX_EPPM_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>

#include "reconx_utils.hxx"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define MIN3(a, b, c) (MIN(a, MIN(b, c)))
#define MIN4(a, b, c, d) (MIN(a, MIN(b, MIN(c, d))))
#define MAX5(a, b, c, d, e) (MAX(a, MAX(b, MAX(c, MAX(d, e)))))

namespace ReconX {

using namespace std;
using namespace Arith;
using namespace Loop;

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
approx_at_cell_interface(const array<const CCTK_REAL, 5> &gf,
                         const CCTK_INT interface) {
  const CCTK_INT Im = interface;
  const CCTK_INT I = interface + 1;
  const CCTK_INT Ip = interface + 2;
  const CCTK_INT Ipp = interface + 3;

  return 7.0 / 12.0 * (gf[I] + gf[Ip]) - (1.0 / 12.0) * (gf[Im] + gf[Ipp]);
}

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
limit(const array<const CCTK_REAL, 5> &gf, const CCTK_REAL gf_rc,
      const CCTK_INT interface, const CCTK_REAL C,
      const bool &keep_var_positive) {
  const CCTK_INT Im = interface;
  const CCTK_INT I = interface + 1;
  const CCTK_INT Ip = interface + 2;
  const CCTK_INT Ipp = interface + 3;

  if ((MIN(gf[I], gf[Ip]) <= gf_rc) && (gf_rc <= MAX(gf[I], gf[Ip]))) {
    return gf_rc;
  } else {
    const CCTK_REAL D2a = 3.0 * ((gf[I] + gf[Ip]) - 2.0 * gf_rc);
    const CCTK_REAL D2aL = ((gf[Im] + gf[Ip]) - 2.0 * gf[I]);
    const CCTK_REAL D2aR = ((gf[I] + gf[Ipp]) - 2.0 * gf[Ip]);
    const CCTK_REAL D2aLim = copysign(1.0, D2a) *
                             MIN3(C * fabs(D2aL), C * fabs(D2aR), fabs(D2a)) *
                             1.0 / 3.0;
    const CCTK_REAL gf_avg = 0.5 * (gf[I] + gf[Ip]);
    if (D2a * D2aR < 0 || D2a * D2aL < 0 ||
        (keep_var_positive && D2aLim > gf_avg))
      return gf_avg;
    else
      return gf_avg - D2aLim;
  }
}

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST void
monotonize(const array<const CCTK_REAL, 5> &gf, CCTK_REAL &rc_minus,
           CCTK_REAL &rc_plus, const CCTK_REAL C) {
  const CCTK_INT Imm = 0;
  const CCTK_INT Im = 1;
  const CCTK_INT I = 2;
  const CCTK_INT Ip = 3;
  const CCTK_INT Ipp = 4;

  CCTK_REAL D2aLim = 0;
  CCTK_REAL rhi = 0;
  const CCTK_REAL daplus = rc_plus - gf[I];
  const CCTK_REAL daminus = gf[I] - rc_minus;
  if (daplus * daminus <= 0 || (gf[Imm] - gf[I]) * (gf[I] - gf[Ipp]) <= 0) {
    const CCTK_REAL D2a = -(12.0 * gf[I] - 6.0 * (rc_minus + rc_plus));
    const CCTK_REAL D2aC = (gf[Im] + gf[Ip]) - 2.0 * gf[I];
    const CCTK_REAL D2aL = (gf[Imm] + gf[I]) - 2.0 * gf[Im];
    const CCTK_REAL D2aR = (gf[I] + gf[Ipp]) - 2.0 * gf[Ip];
    if (copysign(1.0, D2a) == copysign(1.0, D2aC) &&
        copysign(1.0, D2a) == copysign(1.0, D2aL) &&
        copysign(1.0, D2a) == copysign(1.0, D2aR))
      D2aLim = copysign(1.0, D2a) *
               MIN4(C * fabs(D2aL), C * fabs(D2aR), C * fabs(D2aC), fabs(D2a));
    if (!(fabs(D2a) <= 1e-12 * MAX5(fabs(gf[Imm]), fabs(gf[Im]), fabs(gf[I]),
                                    fabs(gf[Ip]), fabs(gf[Ipp]))))
      rhi = D2aLim / D2a;
    if (!(rhi >= 1.0 - 1e-12)) {
      if (daplus * daminus < 0) {
        rc_plus = gf[I] + daplus * rhi;
        rc_minus = gf[I] - daminus * rhi;
      } else if (fabs(daminus) >= 2.0 * fabs(daplus)) {
        rc_minus = gf[I] - (2.0 * (1.0 - rhi) * daplus + rhi * daminus);
      } else if (fabs(daplus) >= 2.0 * fabs(daminus)) {
        rc_plus = gf[I] + (2.0 * (1.0 - rhi) * daminus + rhi * daplus);
      }
    }
  } else {
    if (fabs(daplus) >= 2.0 * fabs(daminus))
      rc_plus = gf[I] + 2.0 * daminus;
    else if (fabs(daminus) >= 2.0 * fabs(daplus))
      rc_minus = gf[I] - 2.0 * daplus;
  }

  return;
}

/* ePPM reconstruction scheme. (see Reisswig et al. 2013) based on McCorquodale
   & Colella (2011) */
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<CCTK_REAL, 2>
eppm(const GF3D2<const CCTK_REAL> &gf_var,
     const array<const vect<int, dim>, 5> &cells, const bool &keep_var_positive,
     const GF3D2<const CCTK_REAL> &gf_press,
     const GF3D2<const CCTK_REAL> &gf_vel_dir,
     const reconstruct_params_t &reconstruct_params) {
  // Unpack all cells in the stencil
  const auto &Imm = cells[0];
  const auto &Im = cells[1];
  const auto &I = cells[2];
  const auto &Ip = cells[3];
  const auto &Ipp = cells[4];

  const array<const CCTK_REAL, 5> gf_stencil{gf_var(Imm), gf_var(Im), gf_var(I),
                                             gf_var(Ip), gf_var(Ipp)};

  // Unpack all PPM parameters
  const CCTK_REAL &ppm_eps = reconstruct_params.ppm_eps;
  const CCTK_REAL &ppm_small = reconstruct_params.ppm_small;
  const CCTK_REAL &ppm_omega1 = reconstruct_params.ppm_omega1;
  const CCTK_REAL &ppm_omega2 = reconstruct_params.ppm_omega2;
  const CCTK_REAL &enhanced_ppm_C2 = reconstruct_params.enhanced_ppm_C2;

  const int iminus = 0, iplus = 1;

  /* approx at cell interface */
  CCTK_REAL rc_minus = approx_at_cell_interface(gf_stencil, iminus);
  CCTK_REAL rc_plus = approx_at_cell_interface(gf_stencil, iplus);

  /* limit */
  rc_minus =
      limit(gf_stencil, rc_minus, iminus, enhanced_ppm_C2, keep_var_positive);
  rc_plus =
      limit(gf_stencil, rc_plus, iplus, enhanced_ppm_C2, keep_var_positive);

  /* monotonize */
  monotonize(gf_stencil, rc_minus, rc_plus, enhanced_ppm_C2);

  /* apply flattening */
  const CCTK_REAL press_Im = gf_press(Im);
  const CCTK_REAL press_Ip = gf_press(Ip);
  const CCTK_REAL dpress_I = press_Ip - press_Im;
  const CCTK_REAL dpress2 = gf_press(Ipp) - gf_press(Imm);
  const CCTK_REAL w_I = ((fabs(dpress_I) > ppm_eps * MIN(press_Im, press_Ip)) &&
                         (gf_vel_dir(Im) > gf_vel_dir(Ip)))
                            ? 1.0
                            : 0.0;
  const CCTK_REAL ftilde_I =
      (fabs(dpress2) < ppm_small)
          ? 1.0
          : MAX(0.0, 1.0 - w_I * MAX(0.0, ppm_omega2 * (dpress_I / dpress2 -
                                                        ppm_omega1)));
  const CCTK_REAL one_minus_ftildeI_gfI = (1 - ftilde_I) * gf_var(I);
  rc_minus = ftilde_I * rc_minus + one_minus_ftildeI_gfI;
  rc_plus = ftilde_I * rc_plus + one_minus_ftildeI_gfI;

  // Return the lower and upper reconstructed states in cell I
  return array<CCTK_REAL, 2>{rc_minus, rc_plus};
}

} // namespace ReconX

#endif // RECONX_EPPM_HXX
