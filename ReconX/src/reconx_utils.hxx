#ifndef RECONX_UTILS_HXX
#define RECONX_UTILS_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <vec.hxx>

namespace ReconX {
using namespace std;
using namespace Arith;

// Struct used to pass parameters to the reconstruction routine
typedef struct {
  // PPM parameters
  bool ppm_shock_detection, ppm_zone_flattening;
  CCTK_REAL poly_k, poly_gamma;
  CCTK_REAL ppm_eta1, ppm_eta2;
  CCTK_REAL ppm_eps, ppm_eps_shock, ppm_small;
  CCTK_REAL ppm_omega1, ppm_omega2;
  CCTK_REAL enhanced_ppm_C2;
  // WENOZ parameters
  CCTK_REAL weno_eps;
  // MP5 parameters
  CCTK_REAL mp5_alpha;
} reconstruct_params_t;

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

} // namespace ReconX

#endif // RECONX_UTILS_HXX

