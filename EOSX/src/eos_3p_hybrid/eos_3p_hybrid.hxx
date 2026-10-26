#ifndef EOS_3P_HYBRID_HXX
#define EOS_3P_HYBRID_HXX

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>

#include "../eos_3p.hxx"
#include "../eos_1p.hxx"

using namespace std;

// Remember to define umass_
namespace EOSX {

template <class ColdEOS> class eos_3p_hybrid : public eos_3p {
public:
  CCTK_REAL gamma;            // dummy
  CCTK_REAL gamma_th, gm1_th; // thermal gamma and gamma-1
  CCTK_REAL temp_over_eps_th; // ratio between the temperature and the thermal part of the internal energy
  range rgeps;                // eps range parameters (allows eps_min sentinel)
  ColdEOS *eos_c;             // cold EOS (polytrope or piecewise-polytrope)

  // constructor
  CCTK_HOST
  CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_3p_hybrid(
      ColdEOS *eos_c_, CCTK_REAL gamma_th_, CCTK_REAL umass_, const range &rgeps_,
      const range &rgrho_, const range &rgye_)
      : rgeps(rgeps_), eos_c(eos_c_), gamma_th(gamma_th_) {

    if (gamma_th <= 1.0 || gamma_th > 2.0) {
      assert(false);
    }

    gm1_th = gamma_th - 1.0;

    // Ensure subluminal soundspeed and P < E for thermal part (same as ideal gas logic)
    if (gamma_th > 2) {
      // unreachable due to assert above, but keep the intent
      rgeps.max =
          std::min(rgeps.max, CCTK_REAL(1.0) / (gamma_th * (gamma_th - 2)));
    }

    set_range_rho(rgrho_);
    set_range_ye(rgye_);
    temp_over_eps_th = gm1_th * umass_;
    //Strictly speaking should probably convert eps range to eps_th range?
    set_range_temp(range(temp_over_eps_th * rgeps.min, temp_over_eps_th * rgeps.max)); // TODO: Implement temperature;
  }

  // cold helpers
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_cold(const CCTK_REAL rho) const {
    const CCTK_REAL gm1 = eos_c->gm1_from_rho(rho);
    return eos_c->sed_from_gm1(gm1);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  p_cold(const CCTK_REAL rho) const {
    const CCTK_REAL gm1 = eos_c->gm1_from_rho(rho);
    return eos_c->p_from_gm1(gm1);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  hm1_cold(const CCTK_REAL rho) const {
    const CCTK_REAL gm1 = eos_c->gm1_from_rho(rho);
    return eos_c->hm1_from_gm1(gm1);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  cs2_cold(const CCTK_REAL rho) const {
    const CCTK_REAL gm1 = eos_c->gm1_from_rho(rho);
    return eos_c->csnd2_from_gm1(gm1);
  }

  // P(rho, eps) = P_cold(rho) + (Gamma_th-1)*rho*(eps-eps_cold(rho))
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                        const CCTK_REAL ye) const {
    const CCTK_REAL p_c = p_cold(rho);
    const CCTK_REAL eps_c = eps_cold(rho);
    const CCTK_REAL p_th = gm1_th * rho * (eps - eps_c);
    return p_c + p_th;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_rho_press_ye(const CCTK_REAL rho, const CCTK_REAL press,
                        const CCTK_REAL ye) const {
    const CCTK_REAL p_c = p_cold(rho);
    const CCTK_REAL eps_c = eps_cold(rho);
    const CCTK_REAL p_th = press - p_c;
    const CCTK_REAL eps_th = p_th / (gm1_th * rho);
    return eps_c + eps_th;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                       const CCTK_REAL ye) const {
    const CCTK_REAL cs2_c = cs2_cold(rho);
    const CCTK_REAL eps_c = eps_cold(rho);
    const CCTK_REAL h_c = CCTK_REAL(1.0) + hm1_cold(rho);
    const CCTK_REAL eps_th = eps - eps_c;
    const CCTK_REAL h_th = gamma_th * eps_th;
    const CCTK_REAL w = h_th / (h_c + h_th);
    const CCTK_REAL cs2 = (CCTK_REAL(1.0) - w) * cs2_c + w * gm1_th;
    return sqrt(cs2);
  }

  // Temperature not implemented
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                        const CCTK_REAL ye) const {
    CCTK_REAL eps = eps_from_rho_temp_ye(rho, temp, ye);
    return csnd_from_rho_eps_ye(rho, eps, ye);
    /*printf("eos_3p_hybrid: sound speed not implemented\n");
    assert(false);
    return CCTK_REAL(0.0);*/
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  temp_from_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                       const CCTK_REAL ye) const {
    const CCTK_REAL eps_c = eps_cold(rho);  
    const CCTK_REAL eps_th = eps - eps_c;
    return temp_over_eps_th * eps_th;
    /*printf("eos_3p_hybrid: temperature not implemented\n");
    assert(false);
    return CCTK_REAL(0.0);*/
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  temp_from_rho_press_ye(const CCTK_REAL rho, CCTK_REAL &press,
                             const CCTK_REAL ye) const {
    printf("eos_3p_hybrid: temperature not implemented");
    return 0.0;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  rho_from_press_temp_ye(CCTK_REAL &press, const CCTK_REAL temp,
                             const CCTK_REAL ye) const {
    printf("eos_3p_hybrid: temperature not implemented");
    return 0.0;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  press_derivs_from_rho_eps_ye(CCTK_REAL &press, CCTK_REAL &dpdrho,
                               CCTK_REAL &dpdeps, const CCTK_REAL rho,
                               const CCTK_REAL eps,
                               const CCTK_REAL ye) const {
    const CCTK_REAL p_c = p_cold(rho);
    const CCTK_REAL cs2_c = cs2_cold(rho);
    const CCTK_REAL eps_c = eps_cold(rho);
    const CCTK_REAL h_c = CCTK_REAL(1.0) + hm1_cold(rho);
    const CCTK_REAL eps_th = eps - eps_c;
    const CCTK_REAL p_th = gm1_th * rho * eps_th;

    press = p_c + p_th;
    dpdrho = h_c * cs2_c + gm1_th * (eps_th - p_c / rho);
    dpdeps = gm1_th * rho;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  entropy_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                           const CCTK_REAL ye) const {
    const CCTK_REAL eps = eps_from_rho_temp_ye(rho, temp, ye);
    return entropy_from_rho_eps_ye(rho, eps, ye);
    /*printf("EOS: entropy from temperature not implemented for eos_3p_hybrid.\n");
    assert(false);
    return CCTK_REAL(0.0);*/
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  entropy_from_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps,
                          const CCTK_REAL ye) const {
    const CCTK_REAL eps_th = eps - eps_cold(rho);
    return log(eps_th * pow(rho, -gm1_th));
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                       const CCTK_REAL ye) const {
    const CCTK_REAL eps_c = eps_cold(rho);
    return temp / temp_over_eps_th + eps_c;
    /* printf("EOS: eps from temperature not implemented for eos_3p_hybrid.\n");
    assert(false);
    return CCTK_REAL(0.0); */
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                         const CCTK_REAL ye) const {
    CCTK_REAL eps = eps_from_rho_temp_ye(rho, temp, ye);    
    return press_from_rho_eps_ye(rho, eps, ye);
    /*printf("EOS: press from temperature not implemented for eos_3p_hybrid.\n");
    assert(false);
    return CCTK_REAL(0.0);*/
  }

  // kappa is the evolved entropy-like quantity.
  // For the *thermal* part: kappa = p_th * rho^(-gamma_th).
  // Hence p_th = kappa * rho^(gamma_th)
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_rho_kappa_ye(const CCTK_REAL rho, const CCTK_REAL kappa,
                          const CCTK_REAL ye) const {
    const CCTK_REAL p_th = kappa * pow(rho, gamma_th);
    return p_cold(rho) + p_th;
  }

  // eps_th = p_th / ((gamma_th-1)*rho) = kappa*rho^(gamma_th-1)/(gamma_th-1)
  // eps = eps_cold + eps_th
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_rho_kappa_ye(const CCTK_REAL rho, const CCTK_REAL kappa,
                        const CCTK_REAL ye) const {
    const CCTK_REAL eps_th =
        kappa * pow(rho, gamma_th - CCTK_REAL(1.0)) / gm1_th;
    return eps_cold(rho) + eps_th;
  }

  // kappa = p_th * rho^(-gamma_th) = (gm1_th*rho*eps_th)*rho^(-gamma_th)
  //       = gm1_th * eps_th * rho^(1-gamma_th)
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  kappa_from_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                        const CCTK_REAL ye) const {
    const CCTK_REAL eps_th = eps - eps_cold(rho);
    return gm1_th * eps_th * pow(rho, CCTK_REAL(1.0) - gamma_th);
  }

  //   eps_min >= 0  -> constant eps_min
  //   eps_min == -1 -> eps_min = eps_cold(rho)
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_3p::range
  range_eps_from_rho_ye(const CCTK_REAL rho, const CCTK_REAL ye) const {
    if (rgeps.min >= CCTK_REAL(0.0)) {
      return range(rgeps.min, rgeps.max);
    } else if (rgeps.min == CCTK_REAL(-1.0)) {
      return range(eps_cold(rho), rgeps.max);
    } else {
      assert(false);
      return range(CCTK_REAL(0.0), rgeps.max);
    }
  }
};

} // namespace EOSX

#endif

