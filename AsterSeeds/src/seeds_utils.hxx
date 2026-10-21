#ifndef SEEDS_UTILS_HXX
#define SEEDS_UTILS_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include "aster_utils.hxx"

namespace AsterSeeds {
using namespace std;
using namespace Loop;
using namespace AsterUtils;

template <typename EOSType>
CCTK_HOST CCTK_DEVICE inline bool
set_beta_floor(const EOSType *eos, const CCTK_REAL press_min,
                const CCTK_REAL temp, const CCTK_REAL Ye, CCTK_REAL &rho,
                CCTK_REAL &eps, CCTK_REAL &press, CCTK_REAL &entropy) {
  if (!std::isfinite(press_min) || press_min < 0.0 ||
      !std::isfinite(temp) || temp < eos->rgtemp.min || temp > eos->rgtemp.max ||
      !std::isfinite(Ye) || Ye < eos->rgye.min || Ye > eos->rgye.max)
    return false;
  CCTK_REAL pressL = press_min;
  CCTK_REAL rhoL = eos->rho_from_press_temp_ye(pressL, temp, Ye);
  // The inverse updates pressL when the target is outside the EOS domain.
  const CCTK_REAL tol = 64.0 * std::numeric_limits<CCTK_REAL>::epsilon();
  if (!std::isfinite(rhoL) || rhoL < (1.0 - tol) * eos->rgrho.min ||
      rhoL > (1.0 + tol) * eos->rgrho.max ||
      !std::isfinite(pressL) || pressL < (1.0 - tol) * press_min)
    return false;
  rhoL = std::clamp(rhoL, eos->rgrho.min, eos->rgrho.max);
  const CCTK_REAL epsL = eos->eps_from_rho_temp_ye(rhoL, temp, Ye);
  const CCTK_REAL entL = eos->kappa_from_rho_temp_ye(rhoL, temp, Ye);
  pressL = eos->press_from_rho_temp_ye(rhoL, temp, Ye);
  if (!std::isfinite(epsL) || !std::isfinite(entL) ||
      !std::isfinite(pressL) || pressL < 0.0)
    return false;
  rho = rhoL;
  eps = epsL;
  press = pressL;
  entropy = entL;
  return true;
}

} // namespace AsterSeeds

#endif // SEEDS_UTILS_HXX
