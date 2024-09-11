#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

#include "aster_utils.hxx"

struct metric {
  CCTK_REAL gxx, gxy, gxz, gyy, gyz, gzz;
};

namespace AsterX {
using namespace Loop;
using namespace Arith;
using namespace AsterUtils;

template <int dir> void ComputeStaggeredB(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputedBstagFromA;
  DECLARE_CCTK_PARAMETERS;

  const std::array dx{CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1),
                      CCTK_DELTA_SPACE(2)};
  const std::array idx{1 / dx[0], 1 / dx[1], 1 / dx[2]};

  static_assert(dir >= 0 && dir < 3, "");

  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};
  grid.loop_int_device<face_centred[0], face_centred[1], face_centred[2]>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" cell indices
        const auto ipjk = p.I + p.DI[0];
        const auto ijpk = p.I + p.DI[1];
        const auto ijkp = p.I + p.DI[2];

        if (dir == 0) {
          /* dBx is curl(A) at (i-1/2,j,k) */
          dBx_stag(p.I) = idx[1] * (Avec_z(ijpk) - Avec_z(p.I)) -
                          idx[2] * (Avec_y(ijkp) - Avec_y(p.I));
        } else if (dir == 1) {
          /* dBy is curl(A) at (i,j-1/2,k) */
          dBy_stag(p.I) = idx[2] * (Avec_x(ijkp) - Avec_x(p.I)) -
                          idx[0] * (Avec_z(ipjk) - Avec_z(p.I));
        } else if (dir == 2) {
          /* dBz is curl(A) at (i,j,z-1/2) */
          dBz_stag(p.I) = idx[0] * (Avec_y(ipjk) - Avec_y(p.I)) -
                          idx[1] * (Avec_x(ijpk) - Avec_x(p.I));
        }

        // TODO: need to implement copy conditions?
      });
}

extern "C" void AsterX_ComputedBstagFromA(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputedBstagFromA;
  DECLARE_CCTK_PARAMETERS;

  ComputeStaggeredB<0>(cctkGH);
  ComputeStaggeredB<1>(cctkGH);
  ComputeStaggeredB<2>(cctkGH);
}

extern "C" void AsterX_ComputedBFromdBstag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputedBFromdBstag;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" cell indices
        const auto ipjk = p.I + p.DI[0];
        const auto ijpk = p.I + p.DI[1];
        const auto ijkp = p.I + p.DI[2];

        /* Second order interpolation of staggered B components to cell center
         */
        dBx(p.I) = 0.5 * (dBx_stag(p.I) + dBx_stag(ipjk));
        dBy(p.I) = 0.5 * (dBy_stag(p.I) + dBy_stag(ijpk));
        dBz(p.I) = 0.5 * (dBz_stag(p.I) + dBz_stag(ijkp));
      });
}

extern "C" void AsterX_ComputeBFromdB(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputeBFromdB;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* Interpolate metric terms from vertices to center */
        metric g;
        g.gxx = calc_avg_v2c(gxx, p);
        g.gxy = calc_avg_v2c(gxy, p);
        g.gxz = calc_avg_v2c(gxz, p);
        g.gyy = calc_avg_v2c(gyy, p);
        g.gyz = calc_avg_v2c(gyz, p);
        g.gzz = calc_avg_v2c(gzz, p);

        /* Determinant of spatial metric */
        const smat<CCTK_REAL, 3> gmat{g.gxx, g.gxy, g.gxz, g.gyy, g.gyz, g.gzz};
        const CCTK_REAL sqrt_detg = sqrt(calc_det(gmat));

        /* Second order interpolation of staggered B components to cell center
         */
        Bvecx(p.I) = dBx(p.I) / sqrt_detg;
        Bvecy(p.I) = dBy(p.I) / sqrt_detg;
        Bvecz(p.I) = dBz(p.I) / sqrt_detg;
      });
}

} // namespace AsterX
