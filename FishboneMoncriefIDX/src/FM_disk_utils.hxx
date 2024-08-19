#ifndef FM_UTILS_HXX
#define FM_UTILS_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace FM_Utils {
using namespace std;
using namespace Loop;

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow2(T x) {
  return x * x;
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

} // namespace FM_Utils

#endif // FM_UTILS_HXX
