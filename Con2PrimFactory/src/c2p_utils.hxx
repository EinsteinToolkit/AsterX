#ifndef C2PUTILS_HXX
#define C2PUTILS_HXX

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
#include <loop_device.hxx>

namespace Con2PrimFactory {
using namespace std;
using namespace Loop;
using namespace Arith;

enum class ROOTSTAT {
  SUCCESS,
  NOT_CONVERGED,
  NOT_BRACKETED
};

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
  return sum_symm<D>([&](int i, int j)
                         ARITH_INLINE { return g(i, j) * v1(i) * v2(j); });
}

template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<T, 3>
calc_cross_product(const vec<T, 3> &B, const vec<T, 3> &v) {
  return ([&](int i) ARITH_INLINE {
    const int j = (i + 1) % 3, k = (i + 2) % 3;
    return B(j) * v(k) - B(k) * v(j);
  });
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

template <typename T, int F>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<vec<T, F>, 3>
calc_cross_product(const vec<vec<T, F>, 3> &B_rc,
                   const vec<vec<T, F>, 3> &v_rc) {
  return ([&](int i) ARITH_INLINE {
    const int j = (i + 1) % 3, k = (i + 2) % 3;
    return vec<T, F>([&](int f) ARITH_INLINE {
      return B_rc(j)(f) * v_rc(k)(f) - B_rc(k)(f) * v_rc(j)(f);
    });
  });
}

template <typename T, int F>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<vec<T, F>, 3>
calc_cross_product(const vec<T, 3> &B, const vec<vec<T, F>, 3> &v_rc) {
  return ([&](int i) ARITH_INLINE {
    const int j = (i + 1) % 3, k = (i + 2) % 3;
    return vec<T, F>([&](int f) ARITH_INLINE {
      return B(j) * v_rc(k)(f) - B(k) * v_rc(j)(f);
    });
  });
}

template <typename T, int F>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<T, F>
pow2(vec<T, F> x) {
  return ([&](int f) { return x(f) * x(f); });
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

template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T pow2(T x) {
  return x * x;
}

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
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd2_v2c(const GF3D2<const T> &gf, const PointDesc &p, int dir) {
  T dgf1, dgf2, dgf3, dgf4;
  int dir1, dir2;
  if (dir == 0) {
    dir1 = 1;
    dir2 = 2;
  } else if (dir == 1) {
    dir1 = 0;
    dir2 = 2;
  } else {
    dir1 = 0;
    dir2 = 1;
  }

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
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd4_v2c(const GF3D2<const T> &gf, const PointDesc &p, int dir) {
  T dgf1, dgf2, dgf3, dgf4;
  T dgf5, dgf6, dgf7, dgf8;
  int dir1, dir2;

  if (dir == 0) {
    dir1 = 1;
    dir2 = 2;
  } else if (dir == 1) {
    dir1 = 0;
    dir2 = 2;
  } else {
    dir1 = 0;
    dir2 = 1;
  }

  dgf1 =
      ((1.0 / 24.0) * gf(p.I + 2 * p.DI[dir]) - (9.0 / 8.0) * gf(p.I + p.DI[dir]) +
       (9.0 / 8.0) * gf(p.I) - (1.0 / 24.0) * gf(p.I - p.DI[dir])) /
      p.DX[dir];
  dgf2 = ((1.0 / 24.0) * gf(p.I + p.DI[dir1] + 2 * p.DI[dir]) -
          (9.0 / 8.0) * gf(p.I + p.DI[dir1] + p.DI[dir]) +
          (9.0 / 8.0) * gf(p.I + p.DI[dir1]) -
          (1.0 / 24.0) * gf(p.I + p.DI[dir1] - p.DI[dir])) /
         p.DX[dir];
  dgf3 = ((1.0 / 24.0) * gf(p.I + p.DI[dir2] + 2 * p.DI[dir]) -
          (9.0 / 8.0) * gf(p.I + p.DI[dir2] + p.DI[dir]) +
          (9.0 / 8.0) * gf(p.I + p.DI[dir2]) -
          (1.0 / 24.0) * gf(p.I + p.DI[dir2] - p.DI[dir])) /
         p.DX[dir];
  dgf4 = ((1.0 / 24.0) * gf(p.I + p.DI[dir1] + p.DI[dir2] + 2 * p.DI[dir]) -
          (9.0 / 8.0) * gf(p.I + p.DI[dir1] + p.DI[dir2] + p.DI[dir]) +
          (9.0 / 8.0) * gf(p.I + p.DI[dir1] + p.DI[dir2]) -
          (1.0 / 24.0) * gf(p.I + p.DI[dir1] + p.DI[dir2] - p.DI[dir])) /
         p.DX[dir];

  dgf5 = ((1.0 / 24.0) * gf(p.I + 2 * p.DI[dir1] + 2 * p.DI[dir2] + 2 * p.DI[dir]) -
          (9.0 / 8.0) * gf(p.I + 2 * p.DI[dir1] + 2 * p.DI[dir2] + p.DI[dir]) +
          (9.0 / 8.0) * gf(p.I + 2 * p.DI[dir1] + 2 * p.DI[dir2]) -
          (1.0 / 24.0) * gf(p.I + 2 * p.DI[dir1] + 2 * p.DI[dir2] - p.DI[dir])) /
         p.DX[dir];

  dgf6 = ((1.0 / 24.0) * gf(p.I + 2 * p.DI[dir1] - p.DI[dir2] + 2 * p.DI[dir]) -
          (9.0 / 8.0) * gf(p.I + 2 * p.DI[dir1] - p.DI[dir2] + p.DI[dir]) +
          (9.0 / 8.0) * gf(p.I + 2 * p.DI[dir1] - p.DI[dir2]) -
          (1.0 / 24.0) * gf(p.I + 2 * p.DI[dir1] - p.DI[dir2] - p.DI[dir])) /
         p.DX[dir];

  dgf7 = ((1.0 / 24.0) * gf(p.I - p.DI[dir1] + 2 * p.DI[dir2] + 2 * p.DI[dir]) -
          (9.0 / 8.0) * gf(p.I - p.DI[dir1] + 2 * p.DI[dir2] + p.DI[dir]) +
          (9.0 / 8.0) * gf(p.I - p.DI[dir1] + 2 * p.DI[dir2]) -
          (1.0 / 24.0) * gf(p.I - p.DI[dir1] + 2 * p.DI[dir2] - p.DI[dir])) /
         p.DX[dir];

  dgf8 = ((1.0 / 24.0) * gf(p.I - p.DI[dir1] - p.DI[dir2] + 2 * p.DI[dir]) -
          (9.0 / 8.0) * gf(p.I - p.DI[dir1] - p.DI[dir2] + p.DI[dir]) +
          (9.0 / 8.0) * gf(p.I - p.DI[dir1] - p.DI[dir2]) -
          (1.0 / 24.0) * gf(p.I - p.DI[dir1] - p.DI[dir2] - p.DI[dir])) /
         p.DX[dir];

  return (27.0 / 14.0) * (dgf1 + dgf2 + dgf3 + dgf4) -
         (10.0 / 7.0) * (dgf5 + dgf6 + dgf7 + dgf8);
}

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

  for (int dk = 0; dk < (dir == 2 ? 1 : 2); ++dk) {
    for (int dj = 0; dj < (dir == 1 ? 1 : 2); ++dj) {
      for (int di = 0; di < (dir == 0 ? 1 : 2); ++di) {
        gf_avg += gf(p.I + p.DI[0] * di + p.DI[1] * dj + p.DI[2] * dk);
      }
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

  for (int dk = 0; dk < (dir == 2 ? 1 : 2); ++dk) {
    for (int dj = 0; dj < (dir == 1 ? 1 : 2); ++dj) {
      for (int di = 0; di < (dir == 0 ? 1 : 2); ++di) {
        gf_avg += gf(p.I + p.DI[0] * di + p.DI[1] * dj + p.DI[2] * dk);
      }
    }
  }
  return gf_avg / 4.0;
}

// Second-order average of cell-centered grid functions to edge center
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_c2e(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  T gf_avg = 0.0;

  for (int dk = 0; dk < (dir == 2 ? 1 : 2); ++dk) {
    for (int dj = 0; dj < (dir == 1 ? 1 : 2); ++dj) {
      for (int di = 0; di < (dir == 0 ? 1 : 2); ++di) {
        gf_avg += gf(p.I - p.DI[0] * di - p.DI[1] * dj - p.DI[2] * dk);
      }
    }
  }
  return gf_avg / 4.0;
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

} // namespace Con2PrimFactory

#endif // #ifndef C2PUTILS_HXX
