#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include "prim2con.hxx"

namespace AsterX {
using namespace Loop;
using namespace std;

extern "C" void AsterX_Prim2Con_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Prim2Con_Initial;
  DECLARE_CCTK_PARAMETERS;

  // Loop over the entire grid (0 to n-1 cells in each direction)
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Interpolate metric terms from vertices to center
        const smat<CCTK_REAL, 3> g{calc_avg_v2c(gxx, p), calc_avg_v2c(gxy, p),
                                   calc_avg_v2c(gxz, p), calc_avg_v2c(gyy, p),
                                   calc_avg_v2c(gyz, p), calc_avg_v2c(gzz, p)};

        prim pv;
        pv.rho = rho(p.I);
        pv.vel(0) = velx(p.I);
        pv.vel(1) = vely(p.I);
        pv.vel(2) = velz(p.I);
        pv.eps = eps(p.I);
        pv.press = press(p.I);
        pv.entropy = entropy(p.I);
        pv.Ye = Ye(p.I);
        pv.Bvec(0) = Bvecx(p.I);
        pv.Bvec(1) = Bvecy(p.I);
        pv.Bvec(2) = Bvecz(p.I);

        cons cv;
        prim2con(g, pv, cv);

        dens(p.I) = cv.dens;
        momx(p.I) = cv.mom(0);
        momy(p.I) = cv.mom(1);
        momz(p.I) = cv.mom(2);
        tau(p.I) = cv.tau;
        DYe(p.I) = cv.DYe;
        dBx(p.I) = cv.dBvec(0);
        dBy(p.I) = cv.dBvec(1);
        dBz(p.I) = cv.dBvec(2);
        DEnt(p.I) = cv.DEnt;
      });
}

extern "C" void AsterX_PsiZero_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_PsiZero_Initial;
  DECLARE_CCTK_PARAMETERS;

  /* Initilaize Psi to 0.0 */
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Psi(p.I) = 0.0; });
}

extern "C" void AsterX_FluxAuxZero_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_FluxAuxZero_Initial;
  DECLARE_CCTK_PARAMETERS;

  /* Auxiliary variables for rhs of Avec and Psi */
  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      Fx(p.I) = 0.0;
                                      Fy(p.I) = 0.0;
                                      Fz(p.I) = 0.0;
                                      Fbetax(p.I) = 0.0;
                                      Fbetay(p.I) = 0.0;
                                      Fbetaz(p.I) = 0.0;
                                      G(p.I) = 0.0;
                                    });

#if 0
  /* Initilaize Flux to 0.0 */
  grid.loop_all_device<0, 1, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      fxdens(p.I) = 0.0;
                                      fxmomx(p.I) = 0.0;
                                      fxmomy(p.I) = 0.0;
                                      fxmomz(p.I) = 0.0;
                                      fxtau(p.I) = 0.0;
                                      fxDEnt(p.I) = 0.0;
                                      fxDYe(p.I) = 0.0;
                                      fxBx(p.I) = 0.0;
                                      fxBy(p.I) = 0.0;
                                      fxBz(p.I) = 0.0;
                                      vtilde_y_xface(p.I) = 0.0;
                                      vtilde_z_xface(p.I) = 0.0;
                                      amax_xface(p.I) = 0.0;
                                      amin_xface(p.I) = 0.0;
                                    });

  grid.loop_all_device<1, 0, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      fydens(p.I) = 0.0;
                                      fymomx(p.I) = 0.0;
                                      fymomy(p.I) = 0.0;
                                      fymomz(p.I) = 0.0;
                                      fytau(p.I) = 0.0;
                                      fyDEnt(p.I) = 0.0;
                                      fyDYe(p.I) = 0.0;
                                      fyBx(p.I) = 0.0;
                                      fyBy(p.I) = 0.0;
                                      fyBz(p.I) = 0.0;
                                      vtilde_x_yface(p.I) = 0.0;
                                      vtilde_z_yface(p.I) = 0.0;
                                      amax_yface(p.I) = 0.0;
                                      amin_yface(p.I) = 0.0;
                                    });

  grid.loop_all_device<1, 1, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      fzdens(p.I) = 0.0;
                                      fzmomx(p.I) = 0.0;
                                      fzmomy(p.I) = 0.0;
                                      fzmomz(p.I) = 0.0;
                                      fztau(p.I) = 0.0;
                                      fzDEnt(p.I) = 0.0;
                                      fzDYe(p.I) = 0.0;
                                      fzBx(p.I) = 0.0;
                                      fzBy(p.I) = 0.0;
                                      fzBz(p.I) = 0.0;
                                      vtilde_x_zface(p.I) = 0.0;
                                      vtilde_y_zface(p.I) = 0.0;
                                      amax_zface(p.I) = 0.0;
                                      amin_zface(p.I) = 0.0;
                                    });
#endif
}

} // namespace AsterX
