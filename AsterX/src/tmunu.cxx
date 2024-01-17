#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "utils.hxx"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;

extern "C" void AsterX_Tmunu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Tmunu;
  DECLARE_CCTK_PARAMETERS;

  /* There are at least two strategies to update vertex-centered Tmunu
   * 1) Use cell-centered mhd variables to first compute cell-centered Tmunu
   * grid-functions, then interpolate cell-centered Tmunu to vertex-centered
   * Tmunu.
   * 2) First, interpolate all mhd quantities to the vertex and then
   * update vertex-centered Tmunu.
   *
   * NOTE: Here, we adopt the (2) strategy.
   * TODO: Include higher order interpolation
   *
   */

  /* grid functions */
  const vec<GF3D2<const CCTK_REAL>, dim> gf_vels{velx, vely, velz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_Bvecs{Bvecx, Bvecy, Bvecz};

  /* Loop over vertex-centers for the entire grid (0 to n-1 cells in each
   * direction) */
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* Interpolating mhd quantities to vertices */

        const CCTK_REAL rho_avg = calc_avg_c2v_4th(rho, p);
        const vec<CCTK_REAL, 3> vup_avg([&](int i) ARITH_INLINE {
          return calc_avg_c2v_4th(gf_vels(i), p);
        });
        const CCTK_REAL eps_avg = calc_avg_c2v_4th(eps, p);
        const CCTK_REAL press_avg = calc_avg_c2v_4th(press, p);
        const vec<CCTK_REAL, 3> Bup_avg([&](int i) ARITH_INLINE {
          return calc_avg_c2v_4th(gf_Bvecs(i), p);
        });

        const smat<CCTK_REAL, 3> g_low{gxx(p.I), gxy(p.I), gxz(p.I),
                                       gyy(p.I), gyz(p.I), gzz(p.I)};
        const vec<CCTK_REAL, 3> beta_up{betax(p.I), betay(p.I), betaz(p.I)};

        /* Computing vlow */
        const vec<CCTK_REAL, 3> vlow_avg = calc_contraction(g_low, vup_avg);

        /* Computing betalow */
        const vec<CCTK_REAL, 3> beta_low = calc_contraction(g_low, beta_up);

        /* Computing beta_sq */
        const CCTK_REAL beta_sq = calc_contraction(beta_low, beta_up);

        /* Computing Lorentz factor */
        const CCTK_REAL w_lor = calc_wlorentz(vup_avg, vlow_avg);

        /* Computing [ \rho(1+\epsilon) + Pgas ]*W^2 = \rho * h * W^2 */
        const CCTK_REAL rhoenthalpyW2 =
            (rho_avg * (1.0 + eps_avg) + press_avg) * w_lor * w_lor;

        /* Computing lower components of 4-velocity (without the Lorentz factor)
         */
        const CCTK_REAL ut_low =
            -alp(p.I) + calc_contraction(beta_low, vup_avg);
        const vec<CCTK_REAL, 3> ui_low = vlow_avg;

        /* Computing upper components of 4-velocity (without the Lorentz factor)
         */
        // utup = 1/alp(p.I); //not used
        const vec<CCTK_REAL, 3> ui_up = vup_avg - beta_up / alp(p.I);

        /* Computing the upper 4 vector b of the magnetic field */
        const CCTK_REAL bst_up =
            w_lor * (calc_contraction(Bup_avg, vlow_avg)) / alp(p.I);
        const vec<CCTK_REAL, 3> bsi_up =
            (Bup_avg + alp(p.I) * bst_up * w_lor * ui_up) / w_lor;

        /* Computing the lower 4 vector b of the magnetic field */
        const CCTK_REAL bst_low = bst_up * (-alp(p.I) * alp(p.I) + beta_sq) +
                                  calc_contraction(bsi_up, beta_low);
        const vec<CCTK_REAL, 3> bsi_low =
            bst_up * beta_low + calc_contraction(g_low, bsi_up);

        /* Calculating b^2 = b^mu b_mu */
        const CCTK_REAL bs2 =
            bst_up * bst_low + calc_contraction(bsi_up, bsi_low);

        /* Computing lower T_{\mu\nu} */
        const CCTK_REAL t00 = rhoenthalpyW2 * ut_low * ut_low +
                              press_avg * (-alp(p.I) * alp(p.I) + beta_sq) +
                              (w_lor * w_lor * ut_low * ut_low +
                               0.5 * (-alp(p.I) * alp(p.I) + beta_sq)) *
                                  bs2 -
                              bst_low * bst_low;
        const vec<CCTK_REAL, 3> t0 =
            rhoenthalpyW2 * ut_low * ui_low + press_avg * beta_low +
            (w_lor * w_lor * ut_low * ui_low + 0.5 * beta_low) * bs2 -
            bst_low * bsi_low;

        const smat<CCTK_REAL, 3> t([&](int i, int j) ARITH_INLINE {
          return rhoenthalpyW2 * ui_low(i) * ui_low(j) +
                 press_avg * g_low(i, j) +
                 (w_lor * w_lor * ui_low(i) * ui_low(j) + 0.5 * g_low(i, j)) *
                     bs2 -
                 bsi_low(i) * bsi_low(j);
        });

        /* Update the Tmunu grid functions */
        eTtt(p.I) = eTtt(p.I) + t00;
        eTtx(p.I) = eTtx(p.I) + t0(0);
        eTty(p.I) = eTty(p.I) + t0(1);
        eTtz(p.I) = eTtz(p.I) + t0(2);
        eTxx(p.I) = eTxx(p.I) + t(0, 0);
        eTxy(p.I) = eTxy(p.I) + t(0, 1);
        eTxz(p.I) = eTxz(p.I) + t(0, 2);
        eTyy(p.I) = eTyy(p.I) + t(1, 1);
        eTyz(p.I) = eTyz(p.I) + t(1, 2);
        eTzz(p.I) = eTzz(p.I) + t(2, 2);
      }); // end of loop over grid
}

} // namespace AsterX
