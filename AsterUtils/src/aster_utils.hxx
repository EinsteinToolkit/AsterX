#ifndef ASTER_UTILS_HXX
#define ASTER_UTILS_HXX

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

#include "aster_fd.hxx"
#include "aster_interp.hxx"

namespace AsterUtils {
using namespace std;
using namespace Loop;
using namespace Arith;

// For contractions, make sure to use the correct indices of the vectors
// Computes the contraction of smat and vec
template <typename T, int D>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<T, D>
calc_contraction(const smat<T, D> &g, const vec<T, D> &v) {
  return [&](int i) ARITH_INLINE {
    return sum<D>([&](int j) ARITH_INLINE { return g(i, j) * v(j); });
  };
}

template <typename T, int D>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_contraction(const vec<T, D> &v_up, const vec<T, D> &v_dn) {
  return sum<D>([&](int i) ARITH_INLINE { return v_up(i) * v_dn(i); });
}

template <typename T, int D>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_contraction(const smat<T, D> &g_dn, const smat<T, D> &T_up) {
  return sum_symm<D>([&](int i, int j)
                         ARITH_INLINE { return g_dn(i, j) * T_up(i, j); });
}

template <typename T, int D>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_contraction(const smat<T, D> &g, const vec<T, D> &v1,
                 const vec<T, D> &v2) {
  // return calc_contraction(v2, calc_contraction(g, v1));
  return sum<D>([&](int i, int j)
                    ARITH_INLINE { return g(i, j) * v1(i) * v2(j); });
}

// Contraction for rc case: consider both sides of the face
template <typename T, int D, int F>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<vec<T, F>, D>
calc_contraction(const smat<T, D> &g, const vec<vec<T, F>, D> &v_rc) {
  return ([&](int i) ARITH_INLINE {
    return vec<T, F>([&](int f) ARITH_INLINE {
      return sum<D>([&](int j) ARITH_INLINE { return g(i, j) * v_rc(j)(f); });
    });
  });
}

template <typename T, int D, int F>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<T, F>
calc_contraction(const vec<vec<T, F>, D> &vup_rc,
                 const vec<vec<T, F>, D> &vdn_rc) {
  return ([&](int f) ARITH_INLINE {
    return sum<D>([&](int i)
                      ARITH_INLINE { return vup_rc(i)(f) * vdn_rc(i)(f); });
  });
}

// Cross products
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<T, 3>
calc_cross_product(const vec<T, 3> &B, const vec<T, 3> &v) {
  return ([&](int i) ARITH_INLINE {
    const int j = (i == 0) ? 1 : ((i == 1) ? 2 : 0);
    const int k = (i == 0) ? 2 : ((i == 1) ? 0 : 1);
    return B(j) * v(k) - B(k) * v(j);
  });
}

template <typename T, int F>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<vec<T, F>, 3>
calc_cross_product(const vec<vec<T, F>, 3> &B_rc,
                   const vec<vec<T, F>, 3> &v_rc) {
  return ([&](int i) ARITH_INLINE {
    const int j = (i == 0) ? 1 : ((i == 1) ? 2 : 0);
    const int k = (i == 0) ? 2 : ((i == 1) ? 0 : 1);
    return vec<T, F>([&](int f) ARITH_INLINE {
      return B_rc(j)(f) * v_rc(k)(f) - B_rc(k)(f) * v_rc(j)(f);
    });
  });
}

template <typename T, int F>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<vec<T, F>, 3>
calc_cross_product(const vec<T, 3> &B, const vec<vec<T, F>, 3> &v_rc) {
  return ([&](int i) ARITH_INLINE {
    const int j = (i == 0) ? 1 : ((i == 1) ? 2 : 0);
    const int k = (i == 0) ? 2 : ((i == 1) ? 0 : 1);
    return vec<T, F>([&](int f) ARITH_INLINE {
      return B(j) * v_rc(k)(f) - B(k) * v_rc(j)(f);
    });
  });
}

// Computes the norm of vec, measured with smat
template <typename T, int D>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_norm(const vec<T, D> &v, const smat<T, D> &g) {
  return sum_symm<D>([&](int i, int j)
                         ARITH_INLINE { return g(i, j) * v(i) * v(j); });
}

// Computes the transpose of vec<vec>
template <typename T, int D, int F>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<vec<T, D>, F>
calc_transpose(const vec<vec<T, F>, D> &vv) {
  return ([&](int f) ARITH_INLINE {
    return vec<T, D>([&](int d) ARITH_INLINE { return vv(d)(f); });
  });
}

// Computes the Lorentz factor
template <typename T, int D>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_wlorentz(const vec<T, D> &v_up, const vec<T, D> &v_dn) {
  return 1.0 / sqrt(1.0 - calc_contraction(v_up, v_dn));
}

template <typename T, int D>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_wlorentz_zvec(const vec<T, D> &zvec_up, const vec<T, D> &zvec_dn) {
  return sqrt(1.0 + calc_contraction(zvec_up, zvec_dn));
}

template <typename T, int F>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<T, F>
pow2(vec<T, F> x) {
  return ([&](int f) { return x(f) * x(f); });
}

template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<T, 6>
get_neighbors(const GF3D2<T> &gf, const PointDesc &p) {
  return {gf(p.I - p.DI[0]), gf(p.I + p.DI[0]), gf(p.I - p.DI[1]),
          gf(p.I + p.DI[1]), gf(p.I - p.DI[2]), gf(p.I + p.DI[2])};
}

template <typename T, int D>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_neighbors(const vec<T, D> flag, const vec<T, D> u_nbs,
                   const vec<T, D> u_saved_nbs) {
  return sum<D>([&](int i) ARITH_INLINE {
           return (flag(i) * u_nbs(i) + (1.0 - flag(i)) * u_saved_nbs(i));
         }) /
         CCTK_REAL(D);
}

// Higher order corrections
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
higher_order_correction(const GF3D2<T> &gf, const PointDesc &p, int dir,
                        int correction_order) {
  T correction;
  correction = 0.0;
  if (correction_order == 4) {
    correction = (13.0 / 12.0) * gf(p.I) -
                 (1.0 / 24.0) * (gf(p.I - p.DI[dir]) + gf(p.I + p.DI[dir]));
  } else if (correction_order == 6) {
    correction =
        (1067.0 / 960.0) * gf(p.I) -
        (29.0 / 480.0) * (gf(p.I - p.DI[dir]) + gf(p.I + p.DI[dir])) +
        (3.0 / 640.0) * (gf(p.I - 2 * p.DI[dir]) + gf(p.I + 2 * p.DI[dir]));
  } else {
    correction = gf(p.I);
  }
  return correction;
}

} // namespace AsterUtils

#endif // ASTER_UTILS_HXX
