#ifndef EOS_3P_IDEALGAS_HXX
#define EOS_3P_IDEALGAS_HXX

#include <cmath>

#include "../eos_3p.hxx"

namespace EOSX {

class eos_3p_idealgas : public eos_3p {
public:
  CCTK_REAL gamma, gm1, inv_gamma, temp_over_eps;
  range rgeps;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  init(CCTK_REAL gamma_, CCTK_REAL umass_, range &rgeps_, const range &rgrho_,
       const range &rgye_) {
    gamma = gamma_;
    gm1 = gamma_ - 1;
    inv_gamma = 1 / gamma_;
    rgeps = rgeps_;
    if (gamma < 1) {
      printf("EOS_IdealGas: initialized with gamma < 1. \n");
      assert(false);
    }
    if (gamma > 2) { // Ensure subluminal Soundspeed and P < E
      rgeps.max = std::min(rgeps.max, 1 / (gamma * (gamma - 2)));
    }
    // set_range_h(range(1 + gamma * rgeps.min, 1 + gamma * rgeps.max));
    set_range_rho(rgrho_);
    set_range_ye(rgye_);
    temp_over_eps = gm1 * umass_;
    set_range_temp(range(temp_over_eps * rgeps.min, temp_over_eps * rgeps.max));
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_rho_eps_ye(
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      CCTK_REAL &eps,      ///< Specific internal energy \f$ \epsilon \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const {
    return gm1 * rho * eps;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_rho_press_ye(
      const CCTK_REAL rho,   ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL press, ///< Pressure \f$ P \f$
      const CCTK_REAL ye     ///< Electron fraction \f$ Y_e \f$
  ) const {
    return press / (rho * gm1);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_rho_eps_ye(
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      CCTK_REAL &eps,      ///< Specific internal energy \f$ \epsilon \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const {
    return sqrt(gm1 * eps / (eps + inv_gamma));
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                              const CCTK_REAL ye) const {
    CCTK_REAL eps = eps_from_rho_temp_ye(rho, temp, ye);
    return csnd_from_rho_eps_ye(rho, eps, ye);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  temp_from_rho_eps_ye(
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      CCTK_REAL &eps,      ///< Specific internal energy \f$ \epsilon \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const {
    return temp_over_eps * eps;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  temp_from_rho_press_ye(
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      CCTK_REAL &press,    ///< Pressure \f$ P \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const {
    CCTK_REAL eps = eps_from_rho_press_ye(rho, press, ye);
    return temp_over_eps * eps;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  rho_from_press_temp_ye(
      CCTK_REAL &press,     ///< Pressure \f$ P \f$
      const CCTK_REAL temp, ///< Temperature \f$ T \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const{
    return press * temp_over_eps / (temp * gm1);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  press_derivs_from_rho_eps_ye(
      CCTK_REAL &press,  ///< Pressure \f$ P \f$
      CCTK_REAL &dpdrho, ///< Partial derivative \f$ \frac{\partial P}{\partial
                         ///< \rho} \f$
      CCTK_REAL &dpdeps, ///< Partial derivative \f$ \frac{\partial P}{\partial
                         ///< \epsilon} \f$
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL eps, ///< Specific internal energy \f$ \epsilon \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const {
    press = gm1 * rho * eps;
    dpdrho = gm1 * eps;
    dpdeps = gm1 * rho;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  entropy_from_rho_temp_ye(
      const CCTK_REAL rho,  ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL temp, ///< Temperature \f$ T \f$
      const CCTK_REAL ye    ///< Electron fraction \f$ Y_e \f$
  ) const {
    return log(temp * pow(rho, -gm1) / temp_over_eps);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  entropy_from_rho_eps_ye(
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL eps, ///< Specific internal energy \f$ \epsilon \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const {
    return log(eps * pow(rho, -gm1));
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_rho_temp_ye(
      const CCTK_REAL rho,  ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL temp, ///< Temperature \f$ T \f$ in MeV
      const CCTK_REAL ye    ///< Electron fraction \f$ Y_e \f$
  ) const {
    return temp / temp_over_eps;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_rho_temp_ye(
      const CCTK_REAL rho,  ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL temp, ///< Temperature \f$ T \f$ in MeV
      const CCTK_REAL ye    ///< Electron fraction \f$ Y_e \f$
  ) const {
    CCTK_REAL eps = eps_from_rho_temp_ye(rho, temp, ye);
    return press_from_rho_eps_ye(rho, eps, ye);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_rho_kappa_ye(const CCTK_REAL rho,
                                const CCTK_REAL kappa, // p/rho^gamma
                                const CCTK_REAL ye) const {
    return kappa * pow(rho, gamma);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_rho_kappa_ye(const CCTK_REAL rho,
                              const CCTK_REAL kappa, // p/rho^gamma
                              const CCTK_REAL ye) const {
    return kappa * pow(rho, gamma - 1.0) / (gamma - 1.0);
  };

  // Note that kappa implements a generic thermodynamic quantity
  // meant to describe the "evolved" entropy by an evolution/application
  // thorn.
  // The notion of the "evolved" entropy (kappa) might differ from the
  // definition of the actual entropy (entropy_from..., see above) for different
  // EOS, e.g. for the ideal gas EOS we have kappa = p * (rho)^(-gamma), where
  // gamma is the adiabatic index.
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  kappa_from_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                              const CCTK_REAL ye) const {
    return (gamma - 1.0) * eps * pow(rho, 1.0 - gamma);
  };

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  kappa_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                        const CCTK_REAL ye) const {
    // Evolved kappa is P / rho^gamma, not physical entropy.
    CCTK_REAL eps = eps_from_rho_temp_ye(rho, temp, ye);
    return kappa_from_rho_eps_ye(rho, eps, ye);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline range
  range_eps_from_rho_ye(
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const {
    return rgeps;
  }

  // Note that the function below do not apply for ideal gas EOS
  // They are used for testing nuX with toy problems
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  mu_lepton_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                                   const CCTK_REAL ye) const {
    return 0.0;
  }
};

} // namespace EOSX

#endif
