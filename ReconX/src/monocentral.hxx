#ifndef RECONX_MONOCENTRAL_HXX
#define RECONX_MONOCENTRAL_HXX

#include <cctk.h>

#include "reconx_utils.hxx"

#include <cmath>
#include <array>

namespace ReconX {

using std::array;

template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
monocentral(const T &x, const T &y) {

  using std::fabs;
  using std::min;

  if (sgn(x) != sgn(y))
    return 0;
  else
    return sgn(x) * min(2 * fabs(x), min(2 * fabs(y), fabs(x + y) / 2));
}

template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<T, 2>
monocentral_reconstruct(T q_Imm, T q_Im, T q_Ip, T q_Ipp) {
  // reconstructs values of Im and Ip at the common face between these
  // two cells

  // reconstructed Im on its "plus/right" side
  T var_slope_p{q_Ip - q_Im};
  T var_slope_m{q_Im - q_Imm};
  const T var_m{q_Im + monocentral(var_slope_p, var_slope_m) / 2};

  // reconstructed Ip on its "minus/left" side
  var_slope_p = q_Ipp - q_Ip;
  var_slope_m = q_Ip - q_Im;
  const T var_p{q_Ip - monocentral(var_slope_p, var_slope_m) / 2};

  return {var_m, var_p};
}

} // namespace ReconX

#endif // RECONX_MINMOD_HXX
