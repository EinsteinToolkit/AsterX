/*! \file atmo.hxx
\brief Class definition representing artificial atmosphere.
*/

#ifndef ATMO_HXX
#define ATMO_HXX

#include "prims.hxx"
#include "cons.hxx"

#include "aster_utils.hxx"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <type_traits>

namespace EOSX {
class eos_3p_tabulated3d;
}

namespace Con2PrimFactory {

/// Class representing an artificial atmosphere.
struct atmosphere {
  CCTK_REAL rho_atmo;
  CCTK_REAL eps_atmo;
  CCTK_REAL ye_atmo;
  CCTK_REAL press_atmo;
  CCTK_REAL temp_atmo;
  CCTK_REAL entropy_atmo;
  CCTK_REAL rho_cut;

  CCTK_DEVICE
  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline atmosphere() = default;

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline atmosphere(
      CCTK_REAL rho_, CCTK_REAL eps_, CCTK_REAL Ye_, CCTK_REAL press_,
      CCTK_REAL temp_, CCTK_REAL entropy_, CCTK_REAL rho_cut_)
      : rho_atmo(rho_), eps_atmo(eps_), ye_atmo(Ye_), press_atmo(press_),
        temp_atmo(temp_), entropy_atmo(entropy_), rho_cut(rho_cut_) {}

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline atmosphere &
  operator=(const atmosphere &other) {
    if (this == &other)
      return *this; // Handle self-assignment
    // Copy data members from 'other' to 'this'
    rho_atmo = other.rho_atmo;
    eps_atmo = other.eps_atmo;
    ye_atmo = other.ye_atmo;
    press_atmo = other.press_atmo;
    temp_atmo = other.temp_atmo;
    entropy_atmo = other.entropy_atmo;
    rho_cut = other.rho_cut;
    return *this;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set(prim_vars &pv) const {

    pv.rho = rho_atmo;
    pv.eps = eps_atmo;
    pv.Ye = ye_atmo;
    pv.press = press_atmo;
    pv.temperature = temp_atmo;
    pv.entropy = entropy_atmo;
    pv.vel(0) = 0.0;
    pv.vel(1) = 0.0;
    pv.vel(2) = 0.0;
    pv.w_lor = 1.0;
    pv.E(0) = 0.0;
    pv.E(1) = 0.0;
    pv.E(2) = 0.0;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set(prim_vars &pv, cons_vars &cv, const smat<CCTK_REAL, 3> &g) const {

    set(pv);
    const CCTK_REAL sqrt_detg = sqrt(calc_det(g));
    cv.dens = sqrt_detg * rho_atmo;
    cv.mom(0) = 0.0;
    cv.mom(1) = 0.0;
    cv.mom(2) = 0.0;
    cv.DYe = cv.dens * ye_atmo;
    cv.DEnt = cv.dens * entropy_atmo;
    const vec<CCTK_REAL, 3> &B_up = pv.Bvec;
    const vec<CCTK_REAL, 3> B_low = calc_contraction(g, B_up);
    CCTK_REAL Bsq = calc_contraction(B_up, B_low);
    cv.tau = cv.dens * eps_atmo + 0.5 * sqrt_detg * Bsq;
  }
};

template <typename EOSIDType, typename EOSType>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline atmosphere
make_atmo(const EOSIDType *eos_1p, const EOSType *eos_3p,
          const CCTK_REAL radial_distance, const CCTK_REAL rho_abs_min,
          const CCTK_REAL p_atmo, const CCTK_REAL t_atmo,
          const CCTK_REAL Ye_atmo, const CCTK_REAL r_atmo,
          const CCTK_REAL n_rho_atmo, const CCTK_REAL n_press_atmo,
          const CCTK_REAL n_temp_atmo, const CCTK_REAL atmo_tol,
          const bool thermal_eos_atmo, const bool use_press_atmo,
          const bool Ye_atmo_beq) {
  // Parameters and EOS modes must be validated before calling this helper.
  // Tabulated EOS uses thermal, temperature-primary atmosphere only.
  CCTK_REAL rho_atm =
      radial_distance > r_atmo
          ? rho_abs_min * pow(r_atmo / radial_distance, n_rho_atmo)
          : rho_abs_min;
  rho_atm = std::clamp(rho_atm, eos_3p->rgrho.min, eos_3p->rgrho.max);
  CCTK_REAL Ye_atm =
      std::clamp(Ye_atmo, eos_3p->rgye.min, eos_3p->rgye.max);
  CCTK_REAL eps_atm;
  CCTK_REAL temp_atm;

  if (thermal_eos_atmo && !use_press_atmo) {
    temp_atm = radial_distance > r_atmo
                   ? t_atmo * pow(r_atmo / radial_distance, n_temp_atmo)
                   : t_atmo;
    temp_atm =
        std::clamp(temp_atm, eos_3p->rgtemp.min, eos_3p->rgtemp.max);
    if (Ye_atmo_beq) {
      if constexpr (std::is_same_v<EOSType, EOSX::eos_3p_tabulated3d>)
        Ye_atm = eos_3p->ye_beq_from_rho_temp(rho_atm, temp_atm);
      else
        assert(false); // Parameter validation restricts this to tabulated EOS.
    }
    // Temperature is authoritative; do not independently floor eps.
    eps_atm = eos_3p->eps_from_rho_temp_ye(rho_atm, temp_atm, Ye_atm);
  } else {
    assert(!Ye_atmo_beq);
    if (thermal_eos_atmo) {
      // Pressure-primary atmosphere requires a supported EOS inversion.
      const CCTK_REAL press_atm =
          radial_distance > r_atmo
              ? p_atmo * pow(r_atmo / radial_distance, n_press_atmo)
              : p_atmo;
      eps_atm = eos_3p->eps_from_rho_press_ye(rho_atm, press_atm, Ye_atm);
    } else {
      // Match cold-EOS energy, then use the evolution EOS for closure.
      const CCTK_REAL gm1 = eos_1p->gm1_from_rho(rho_atm);
      eps_atm = eos_1p->sed_from_gm1(gm1);
    }
    const auto rgeps = eos_3p->range_eps_from_rho_ye(rho_atm, Ye_atm);
    eps_atm = std::clamp(eps_atm, rgeps.min, rgeps.max);
    temp_atm = eos_3p->temp_from_rho_eps_ye(rho_atm, eps_atm, Ye_atm);
  }

  // Recompute dependent quantities after bounding the authoritative state.
  const CCTK_REAL press_atm =
      eos_3p->press_from_rho_temp_ye(rho_atm, temp_atm, Ye_atm);
  const CCTK_REAL entropy_atm =
      eos_3p->kappa_from_rho_temp_ye(rho_atm, temp_atm, Ye_atm);
  const CCTK_REAL rho_atmo_cut = rho_atm * (1 + atmo_tol);
  return atmosphere(rho_atm, eps_atm, Ye_atm, press_atm, temp_atm,
                    entropy_atm, rho_atmo_cut);
}

} // namespace Con2PrimFactory
#endif
