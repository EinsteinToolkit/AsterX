#include "prim2con.hxx"

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace AsterX {
using namespace Loop;

extern "C" void AsterX_Prim2Con_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Prim2Con_Initial;
  DECLARE_CCTK_PARAMETERS;

  // Loop over the entire grid (0 to n-1 cells in each direction)
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Interpolate metric terms from vertices to center
        const smat<CCTK_REAL, 3> g{gxx(p.I), gxy(p.I),
                                   gxz(p.I), gyy(p.I),
                                   gyz(p.I), gzz(p.I)};

        prim pv;
        pv.rho = rho(p.I);
        pv.vel(0) = velx(p.I);
        pv.vel(1) = vely(p.I);
        pv.vel(2) = velz(p.I);
        pv.eps = eps(p.I);
        pv.press = press(p.I);
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
        dBx(p.I) = cv.dBvec(0);
        dBy(p.I) = cv.dBvec(1);
        dBz(p.I) = cv.dBvec(2);

        saved_rho(p.I) = pv.rho;
        saved_velx(p.I) = pv.vel(0);
        saved_vely(p.I) = pv.vel(1);
        saved_velz(p.I) = pv.vel(2);
        saved_eps(p.I) = pv.eps;
      });

  /* Initilaize Psi to 0.0 */
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Psi(p.I) = 0.0; });
}

} // namespace AsterX
