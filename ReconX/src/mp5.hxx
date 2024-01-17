#ifndef RECONX_MP5_HXX
#define RECONX_MP5_HXX

#include <cctk.h>

#include "minmod.hxx"

#include <cmath>

namespace ReconX {

// Compute the median of three numbers
template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T median(T &x, T &y,
                                                                   T &z) {
  return x + minmod(y - x, z - x);
}

/**
 * Fifth-order monotonicity preserving (MP5) scheme of Suresh and Huynh
 * 1997. Paper: "Accurate monotonicity-preserving schemes with Runge–Kutta time
 * stepping, J. Comput. Phys. 136 (1997) 83–99."
 */
template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
mp5(T gf_Imm, T gf_Im, T gf_I, T gf_Ip, T gf_Ipp, T mp5_alpha) {

  using std::max;
  using std::min;

  // Eq. (2.1)
  const T ul{(2 * gf_Imm - 13 * gf_Im + 47 * gf_I + 27 * gf_Ip - 3 * gf_Ipp) /
             60.0};

  const T deltam{gf_I - gf_Im};
  const T deltap{gf_Ip - gf_I};

  // Eq. (2.12)
  const T ump{gf_I + minmod(deltap, mp5_alpha * deltam)};

  // Condition Eq. (2.30)
  if ((ul - gf_I) * (ul - ump) <= 1.0e-10) {
    return ul;
  } else {

    // Eq. (2.19)
    const T dm{gf_Imm + gf_I - 2 * gf_Im};
    const T d{gf_Im + gf_Ip - 2 * gf_I};
    const T dp{gf_I + gf_Ipp - 2 * gf_Ip};

    // Eq. (2.27)
    const T dmp{minmod(minmod(4 * d - dp, 4 * dp - d), minmod(d, dp))};
    const T dmm{minmod(minmod(4 * dm - d, 4 * d - dm), minmod(dm, d))};

    // Eq. (2.8)
    const T uul{gf_I + mp5_alpha * deltam};

    // Eq. (2.16)
    const T uav{(gf_I + gf_Ip) / 2.0};

    // Eq. (2.28)
    const T umd{uav - dmp / 2.0};

    // Eq. (2.29)
    const T ulc{gf_I + deltam / 2.0 + 4.0 / 3.0 * dmm};

    // Eq. (2.24 a)
    const T umin{max(min(gf_I, min(gf_Ip, umd)), min(gf_I, min(uul, ulc)))};

    // Eq. (2.24 b)
    const T umax{min(max(gf_I, max(gf_Ip, umd)), max(gf_I, max(uul, ulc)))};

    // Eq. (2.26)
    return median(ul, umin, umax);
  }
  return T{0};
}

template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<T, 2>
mp5_reconstruct(T gf_Immm, T gf_Imm, T gf_Im, T gf_Ip, T gf_Ipp, T gf_Ippp,
                T mp5_alpha) {
  // for the left cell, the plus side has sequence: Immm, Imm, Im, Ip, Ipp
  // for the left cell, the minus side has sequence: Ipp, Ip, Im, Im, Immm
  // here, we need the plus side
  const auto rc_Im{mp5(gf_Immm, gf_Imm, gf_Im, gf_Ip, gf_Ipp, mp5_alpha)};

  // for the right cell, the plus side has sequence: Imm, Im, Ip, Ipp, Ippp
  // for the right cell, the minus side has sequence: Ippp, Ipp, Ip, Im, Imm
  // here, we need the minus side
  const auto rc_Ip{mp5(gf_Ippp, gf_Ipp, gf_Ip, gf_Im, gf_Imm, mp5_alpha)};

  return {rc_Im, rc_Ip};
}

} // namespace ReconX

#endif // RECONX_MP5_HXX
