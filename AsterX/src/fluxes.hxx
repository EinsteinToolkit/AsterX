#ifndef ASTERX_FLUXES_HXX
#define ASTERX_FLUXES_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

namespace AsterX {
using namespace std;
using namespace Arith;

// Lax-Friedrichs solver
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
laxf(vec<vec<CCTK_REAL, 4>, 2> lam, vec<CCTK_REAL, 2> var,
     vec<CCTK_REAL, 2> flux) {
  const CCTK_REAL charmax =
      max({CCTK_REAL(0), fabs(lam(0)(0)), fabs(lam(0)(1)), fabs(lam(0)(2)),
           fabs(lam(0)(3)), fabs(lam(1)(0)), fabs(lam(1)(1)), fabs(lam(1)(2)),
           fabs(lam(1)(3))});

  return 0.5 * ((flux(0) + flux(1)) - charmax * (var(1) - var(0)));
}

// HLLE solver
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
hlle(vec<vec<CCTK_REAL, 4>, 2> lam, vec<CCTK_REAL, 2> var,
     vec<CCTK_REAL, 2> flux) {
  const CCTK_REAL charmax =
      max({CCTK_REAL(0), lam(0)(0), lam(0)(1), lam(0)(2), lam(0)(3), lam(1)(0),
           lam(1)(1), lam(1)(2), lam(1)(3)});

  // Note that charmin is just the minimum, not with the minus sign
  const CCTK_REAL charmin =
      min({CCTK_REAL(0), lam(0)(0), lam(0)(1), lam(0)(2), lam(0)(3), lam(1)(0),
           lam(1)(1), lam(1)(2), lam(1)(3)});

  const CCTK_REAL charpm = charmax - charmin;

  return (charmax * flux(0) - charmin * flux(1) +
          charmax * charmin * (var(1) - var(0))) /
         charpm;
}

} // namespace AsterX

#endif // ASTERX_FLUXES_HXX
