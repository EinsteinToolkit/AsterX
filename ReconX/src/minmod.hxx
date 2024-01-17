#ifndef RECONX_MINMOD_HXX
#define RECONX_MINMOD_HXX

#include <cctk.h>

#include <cmath>
#include <array>

namespace ReconX {

using std::array;

template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T minmod(const T &x,
                                                                   const T &y) {
  using std::fabs;
  using std::signbit;

  if (signbit(x) != signbit(y))
    return T(0);
  if (fabs(x) < fabs(y))
    return x;
  else
    return y;
}

template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<T, 2>
minmod_reconstruct(T q_Imm, T q_Im, T q_Ip, T q_Ipp) {
  // reconstructs values of Im and Ip at the common face between these
  // two cells

  const T var_slope_p{q_Ipp - q_Ip};
  const T var_slope_c{q_Ip - q_Im};
  const T var_slope_m{q_Im - q_Imm};

  // reconstructed Im on its "plus/right" side
  const T var_m{q_Im + minmod(var_slope_c, var_slope_m) / 2};

  // reconstructed Ip on its "minus/left" side
  const T var_p{q_Ip - minmod(var_slope_p, var_slope_c) / 2};

  return {var_m, var_p};
}

} // namespace ReconX

#endif // RECONX_MINMOD_HXX
