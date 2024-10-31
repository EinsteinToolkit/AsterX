#ifndef ASTER_FD_HXX
#define ASTER_FD_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <mat.hxx>
#include <vec.hxx>
#include <sum.hxx>
#include <simd.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>

namespace AsterUtils {
using namespace std;
using namespace Loop;
using namespace Arith;

// FD2: vertex centered input, vertex centered output, oneside stencil
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd2_v2v_oneside(const GF3D2<const T> &gf, const PointDesc &p,
                     const int dir, const int sign) {
  return -sign *
         (gf(p.I + 2 * sign * p.DI[dir]) - 4.0 * gf(p.I + sign * p.DI[dir]) +
          3.0 * gf(p.I)) *
         (0.5 / p.DX[dir]);
}

// FD2: cell centered input, cell centered output
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd2_c2c(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  return (0.5 / p.DX[dir]) * (gf(p.I + p.DI[dir]) - gf(p.I - p.DI[dir]));
}

// FD2: vertex centered input, edge centered output
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd2_v2e(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  return (gf(p.I + p.DI[dir]) - gf(p.I)) / p.DX[dir];
}

// FD2: vertex centered input, cell centered output
template <int FDORDER, typename T>
CCTK_DEVICE CCTK_HOST
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline std::enable_if_t<FDORDER == 2, T>
    calc_fd_v2c(const GF3D2<const T> &gf, const PointDesc &p, int dir) {
  T dgf1, dgf2, dgf3, dgf4;
  const int dir1 = (dir == 0) ? 1 : ((dir == 1) ? 2 : 0);
  const int dir2 = (dir == 0) ? 2 : ((dir == 1) ? 0 : 1);

  dgf1 = (gf(p.I + p.DI[dir]) - gf(p.I)) / p.DX[dir];
  dgf2 = (gf(p.I + p.DI[dir1] + p.DI[dir]) - gf(p.I + p.DI[dir1])) / p.DX[dir];
  dgf3 = (gf(p.I + p.DI[dir2] + p.DI[dir]) - gf(p.I + p.DI[dir2])) / p.DX[dir];
  dgf4 = (gf(p.I + p.DI[dir1] + p.DI[dir2] + p.DI[dir]) -
          gf(p.I + p.DI[dir1] + p.DI[dir2])) /
         p.DX[dir];
  return 0.25 * (dgf1 + dgf2 + dgf3 + dgf4);
}

// FD4: cell centered input, cell centered output
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd4_c2c(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  return (1.0 / (12.0 * p.DX[dir])) *
         (-gf(p.I + 2 * p.DI[dir]) + 8.0 * gf(p.I + p.DI[dir]) -
          8.0 * gf(p.I - p.DI[dir]) + gf(p.I - 2 * p.DI[dir]));
}

// FD4: vertex centered input, cell centered output
// Interpolation from edges to cell centers is 2nd order
template <int FDORDER, typename T>
CCTK_DEVICE CCTK_HOST
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline std::enable_if_t<FDORDER == 4, T>
    calc_fd_v2c(const GF3D2<const T> &gf, const PointDesc &p, int dir) {
  T dgf1, dgf2, dgf3, dgf4;

  const int dir1 = (dir == 0) ? 1 : ((dir == 1) ? 2 : 0);
  const int dir2 = (dir == 0) ? 2 : ((dir == 1) ? 0 : 1);

  dgf1 = ((1.0 / 24.0) * gf(p.I - p.DI[dir]) - (27.0 / 24.0) * gf(p.I) +
          (27.0 / 24.0) * gf(p.I + p.DI[dir]) -
          (1.0 / 24.0) * gf(p.I + 2 * p.DI[dir])) /
         p.DX[dir];

  dgf2 = ((1.0 / 24.0) * gf(p.I + p.DI[dir1] - p.DI[dir]) -
          (27.0 / 24.0) * gf(p.I + p.DI[dir1]) +
          (27.0 / 24.0) * gf(p.I + p.DI[dir1] + p.DI[dir]) -
          (1.0 / 24.0) * gf(p.I + p.DI[dir1] + 2 * p.DI[dir])) /
         p.DX[dir];

  dgf3 = ((1.0 / 24.0) * gf(p.I + p.DI[dir2] - p.DI[dir]) -
          (27.0 / 24.0) * gf(p.I + p.DI[dir2]) +
          (27.0 / 24.0) * gf(p.I + p.DI[dir2] + p.DI[dir]) -
          (1.0 / 24.0) * gf(p.I + p.DI[dir2] + 2 * p.DI[dir])) /
         p.DX[dir];

  dgf4 = ((1.0 / 24.0) * gf(p.I + p.DI[dir1] + p.DI[dir2] - p.DI[dir]) -
          (27.0 / 24.0) * gf(p.I + p.DI[dir1] + p.DI[dir2]) +
          (27.0 / 24.0) * gf(p.I + p.DI[dir1] + p.DI[dir2] + p.DI[dir]) -
          (1.0 / 24.0) * gf(p.I + p.DI[dir1] + p.DI[dir2] + 2 * p.DI[dir])) /
         p.DX[dir];

  return 0.25 * (dgf1 + dgf2 + dgf3 + dgf4);
}

} // namespace AsterUtils

#endif // ASTER_FD_HXX
