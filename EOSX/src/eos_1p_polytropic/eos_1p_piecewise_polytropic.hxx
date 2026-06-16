#ifndef EOS_1P_PIECEWISE_POLYTROPIC_HXX
#define EOS_1P_PIECEWISE_POLYTROPIC_HXX

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "../eos_1p.hxx"

namespace EOSX {

#ifndef EOSX_PWPOLY_MAX_SEGM
#define EOSX_PWPOLY_MAX_SEGM 64
#endif

/// This implements the "cold" part used by the hybrid EOS:
///   P = P_cold(rho),  eps = eps_cold(rho), etc.
///
/// Notation follows EOSX:
///   rho  ... rest-mass density
///   gm1  ... g-1 (= h-1 for isentropic/cold EOS)
///
class eos_1p_piecewise_polytropic : public eos_1p {
public:
  // A single polytropic segment (mirrors whizza::eos_poly_piece)
  struct eos_poly_piece {
    CCTK_REAL rho0;    // segment start density
    CCTK_REAL dsed;    // energy shift to enforce continuity
    CCTK_REAL gamma;   // segment gamma
    CCTK_REAL rho_p;   // polytropic density scale for this segment (rho_p)
    CCTK_REAL n;       // n = 1/(gamma-1)
    CCTK_REAL np1;     // n+1
    CCTK_REAL invn;    // 1/n
    CCTK_REAL gm10;    // gm1 at rho0
    CCTK_REAL p0;      // p at rho0

    eos_poly_piece() = default;

    CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_poly_piece(
        CCTK_REAL rho0_, CCTK_REAL sed0_, CCTK_REAL gamma_, CCTK_REAL rho_p_)
        : rho0(rho0_), gamma(gamma_), rho_p(rho_p_) {
      assert(gamma > 1.0);
      n = CCTK_REAL(1.0) / (gamma - CCTK_REAL(1.0));
      np1 = n + CCTK_REAL(1.0);
      invn = CCTK_REAL(1.0) / n;

      // Continuity shift:
      // sed(rho) = sed_poly(gm1) + dsed, with dsed chosen so sed(rho0)=sed0
      dsed = sed0_ - n * pow(rho0 / rho_p, invn);

      gm10 = gm1_from_rho(rho0);
      p0 = p_from_gm1(gm10);
    }

    // g-1 = h-1 for isentropic EOS; enforce continuity via dsed
    CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
    gm1_from_rho(const CCTK_REAL rho) const {
      // gm1 = (n+1)*(rho/rho_p)^(1/n) + dsed
      return np1 * pow(rho / rho_p, invn) + dsed;
    }

    CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
    gm1_from_p(const CCTK_REAL p) const {
      // gm1 = (n+1)*(p/rho_p)^(1/(n+1)) + dsed
      return np1 * pow(p / rho_p, CCTK_REAL(1.0) / np1) + dsed;
    }

    CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
    sed_from_gm1(const CCTK_REAL gm1) const {
      // sed = (gm1 - dsed)/gamma + dsed
      return (gm1 - dsed) / gamma + dsed;
    }

    CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
    ied_from_gm1(const CCTK_REAL gm1) const {
      // rho_I = rho * sed
      return sed_from_gm1(gm1) * rho_from_gm1(gm1);
    }

    CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
    p_from_gm1(const CCTK_REAL gm1) const {
      // p = rho_p * ((gm1 - dsed)/(n+1))^(n+1)
      return rho_p * pow((gm1 - dsed) / np1, np1);
    }

    CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
    rho_from_gm1(const CCTK_REAL gm1) const {
      // rho = rho_p * ((gm1 - dsed)/(n+1))^n
      return rho_p * pow((gm1 - dsed) / np1, n);
    }

    CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
    hm1_from_gm1(const CCTK_REAL gm1) const {
      return gm1;
    }

    CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
    csnd2_from_gm1(const CCTK_REAL gm1) const {
      // c_s^2 = (gm1 - dsed) / ( n*(gm1 + 1) )
      const CCTK_REAL gm1p = gm1 - dsed;
      return gm1p / (n * (gm1 + CCTK_REAL(1.0)));
    }
  };

  eos_poly_piece segments[EOSX_PWPOLY_MAX_SEGM];
  CCTK_INT nsegm = 0;
  CCTK_REAL rho_max;

  // Initialize from:
  //  - rho_p0: first segment polytropic density scale
  //  - nsegm_: number of segments
  //  - segm_bound: boundary densities (must start at 0), length = nsegm_
  //  - segm_gamma: segment gammas, length = nsegm_
  //  - rho_max_: EOS validity max density
  //
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void init(
      CCTK_REAL rho_p0, const CCTK_INT nsegm_,
      const CCTK_REAL *segm_bound, const CCTK_REAL *segm_gamma,
      CCTK_REAL rho_max_) {

    if (nsegm_ <= 0) {
      printf("eos_1p_piecewise_polytropic: need at least one segment.\n");
      assert(false);
    }
    if (nsegm_ > EOSX_PWPOLY_MAX_SEGM) {
      printf("eos_1p_piecewise_polytropic: nsegm exceeds EOSX_PWPOLY_MAX_SEGM.\n");
      assert(false);
    }
    if (segm_bound[0] != CCTK_REAL(0.0)) {
      printf("eos_1p_piecewise_polytropic: first segment must start at rho=0.\n");
      assert(false);
    }

    nsegm = nsegm_;
    rho_max = rho_max_;

    // First segment anchored at rho=0 with sed0=0
    segments[0] = eos_poly_piece(CCTK_REAL(0.0), CCTK_REAL(0.0),
                                segm_gamma[0], rho_p0);

    // Build subsequent segments enforcing continuity of sed at each boundary
    for (CCTK_INT i = 1; i < nsegm; ++i) {
      const CCTK_REAL rho_b = segm_bound[i];
      const CCTK_REAL gamma_i = segm_gamma[i];

      const eos_poly_piece &lseg = segments[i - 1];
      if (lseg.rho0 >= rho_b) {
        printf("eos_1p_piecewise_polytropic: boundaries not strictly increasing.\n");
        assert(false);
      }

      // sed at boundary from previous segment
      const CCTK_REAL sedc = lseg.sed_from_gm1(lseg.gm1_from_rho(rho_b));

      // Compute rho_p for the new segment so pressure is continuous.
      const CCTK_REAL np = CCTK_REAL(1.0) / (gamma_i - CCTK_REAL(1.0));
      const CCTK_REAL et = np / lseg.n;
      const CCTK_REAL rho_p = pow(lseg.rho_p, et) * pow(rho_b, CCTK_REAL(1.0) - et);

      segments[i] = eos_poly_piece(rho_b, sedc, gamma_i, rho_p);
    }

    // If your eos_1p base has range-setting, this is where it would go.
    // set_ranges(range(0, rho_max));
  }

  // Segment selection utilities
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline const eos_poly_piece &
  segment_for_rho(const CCTK_REAL rho) const {
    // Find last segment with rho0 <= rho
    for (CCTK_INT idx = nsegm - 1; idx >= 0; --idx) {
      if (segments[idx].rho0 <= rho) return segments[idx];
    }
    return segments[0];
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline const eos_poly_piece &
  segment_for_p(const CCTK_REAL p) const {
    // Find last segment with p0 <= p
    for (CCTK_INT idx = nsegm - 1; idx >= 0; --idx) {
      if (segments[idx].p0 <= p) return segments[idx];
    }
    return segments[0];
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline const eos_poly_piece &
  segment_for_gm1(const CCTK_REAL gm1) const {
    // Find last segment with gm10 <= gm1
    for (CCTK_INT idx = nsegm - 1; idx >= 0; --idx) {
      if (segments[idx].gm10 <= gm1) return segments[idx];
    }
    return segments[0];
  }

  // Required cold EOS API (EOSX naming)
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  gm1_from_rho(const CCTK_REAL rho) const {
    return segment_for_rho(rho).gm1_from_rho(rho);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  gm1_from_p(const CCTK_REAL p) const {
    return segment_for_p(p).gm1_from_p(p);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  sed_from_gm1(const CCTK_REAL gm1) const {
    return segment_for_gm1(gm1).sed_from_gm1(gm1);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  ied_from_gm1(const CCTK_REAL gm1) const {
    return segment_for_gm1(gm1).ied_from_gm1(gm1);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  p_from_gm1(const CCTK_REAL gm1) const {
    return segment_for_gm1(gm1).p_from_gm1(gm1);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  rho_from_gm1(const CCTK_REAL gm1) const {
    return segment_for_gm1(gm1).rho_from_gm1(gm1);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  hm1_from_gm1(const CCTK_REAL gm1) const {
    return gm1;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd2_from_gm1(const CCTK_REAL gm1) const {
    return segment_for_gm1(gm1).csnd2_from_gm1(gm1);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  temp_from_gm1(const CCTK_REAL gm1) const {
    const eos_poly_piece &seg = segment_for_gm1(gm1);

    // For each polytropic segment:
    //   gm1 - dsed = (n+1) * P/rho
    // so assuming ideal-gas closure:
    //   theta = mu * P/rho = mu * (gm1 - dsed)/(n+1)
    const CCTK_REAL mu = CCTK_REAL(1.0);

    return mu * (gm1 - seg.dsed) / seg.np1;
  }
  
};

} // namespace EOSX

#endif

