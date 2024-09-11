#ifndef ASTER_INTERP_HXX
#define ASTER_INTERP_HXX

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

namespace AsterUtils {
using namespace std;
using namespace Loop;
using namespace Arith;

// Second-order average of vertex-centered grid functions to cell center
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_v2c(const GF3D2<const T> &gf, const PointDesc &p) {
  T gf_avg = 0.0;

  for (int dk = 0; dk < 2; ++dk) {
    for (int dj = 0; dj < 2; ++dj) {
      for (int di = 0; di < 2; ++di) {
        gf_avg += gf(p.I + p.DI[0] * di + p.DI[1] * dj + p.DI[2] * dk);
      }
    }
  }
  return gf_avg / 8.0;
}

// Second-order average of edge-centered grid functions to vertex-centered
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_e2v(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  T gf_avg = 0.0;

  for (int di = 0; di < 2; ++di) {
    gf_avg += gf(p.I - p.DI[dir] * di);
  }
  return gf_avg / 2.0;
}

// Second-order average of edge-centered grid functions (along dir) to cell
// center
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_e2c(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  T gf_avg = 0.0;
  const int j = (dir == 0) ? 1 : ((dir == 1) ? 2 : 0);
  const int k = (dir == 0) ? 2 : ((dir == 1) ? 0 : 1);
  for (int dk = 0; dk < 2; ++dk) {
    for (int dj = 0; dj < 2; ++dj) {
      gf_avg += gf(p.I + p.DI[j] * dj + p.DI[k] * dk);
    }
  }
  return gf_avg / 4.0;
}

// Second-order average of vertex-centered grid functionsto face
// center (perp to dir)
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_v2f(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  T gf_avg = 0.0;
  const int j = (dir == 0) ? 1 : ((dir == 1) ? 2 : 0);
  const int k = (dir == 0) ? 2 : ((dir == 1) ? 0 : 1);
  for (int dk = 0; dk < 2; ++dk) {
    for (int dj = 0; dj < 2; ++dj) {
      gf_avg += gf(p.I + p.DI[j] * dj + p.DI[k] * dk);
    }
  }
  return gf_avg / 4.0;
}

// Second-order average of cell-centered grid functions to edge center
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_c2e(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  T gf_avg = 0.0;
  const int j = (dir == 0) ? 1 : ((dir == 1) ? 2 : 0);
  const int k = (dir == 0) ? 2 : ((dir == 1) ? 0 : 1);
  for (int dk = 0; dk < 2; ++dk) {
    for (int dj = 0; dj < 2; ++dj) {
      gf_avg += gf(p.I - p.DI[j] * dj - p.DI[k] * dk);
    }
  }
  return gf_avg / 4.0;
}

template <int interp_order, typename T>
CCTK_DEVICE CCTK_HOST
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline std::enable_if_t<interp_order == 2, T>
    calc_avg_c2v(const GF3D2<const T> &gf, const PointDesc &p) {
  T gf_avg = 0.0;

  for (int dk = 0; dk < 2; ++dk) {
    for (int dj = 0; dj < 2; ++dj) {
      for (int di = 0; di < 2; ++di) {
        gf_avg += gf(p.I - p.DI[0] * di - p.DI[1] * dj - p.DI[2] * dk);
      }
    }
  }
  return gf_avg / 8.0;
}

template <int interp_order, typename T>
CCTK_DEVICE CCTK_HOST
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline std::enable_if_t<interp_order == 4, T>
    calc_avg_c2v(const GF3D2<const T> &gf, const PointDesc &p) {
  T gf_avg = 0.0;
  const vect<T, 4> wt = {-1.0 / 16.0, 9.0 / 16.0, 9.0 / 16.0, -1.0 / 16.0};
  for (int i = 0; i < dim; i++) {
    const int j = (i == 0) ? 1 : ((i == 1) ? 2 : 0);
    const int k = (i == 0) ? 2 : ((i == 1) ? 0 : 1);
    // interp in i-dir to get plane x_i = 0
    vect<vect<T, 4>, 4> gf_i = {{0.0, 0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0, 0.0}};
    for (int dk = 0; dk < 4; ++dk) {
      for (int dj = 0; dj < 4; ++dj) {
        for (int di = 0; di < 4; ++di) {
          gf_i[dk][dj] += wt[di] * gf(p.I + p.DI[i] * (di - 2) +
                                      p.DI[j] * (dj - 2) + p.DI[k] * (dk - 2));
        }
      }
    }
    // interp in j-dir to get line x_i = x_j = 0
    vect<T, 4> gf_j = {0.0, 0.0, 0.0, 0.0};
    for (int dk = 0; dk < 4; ++dk) {
      for (int dj = 0; dj < 4; ++dj) {
        gf_j[dk] += wt[dj] * gf_i[dk][dj];
      }
    }
    // interp in k-dir to get point x_i = x_j = x_k = 0
    for (int dk = 0; dk < 4; ++dk) {
      gf_avg += wt[dk] * gf_j[dk];
    }
  }
  return gf_avg / 3.0;
}

} // namespace AsterUtils

#endif // ASTER_INTERP_HXX
