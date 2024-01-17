#ifndef EOS_IDEALGAS_HXX
#define EOS_IDEALGAS_HXX

#include <stdexcept>
#include <algorithm>
#include <cmath>

#include "eos.hxx"
using namespace std;

namespace EOSX {

class eos_idealgas : public eos {
public:
  CCTK_REAL gamma, gm1, temp_over_eps;
  range rgeps;

  // constructor
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_idealgas(
      CCTK_REAL gamma_, CCTK_REAL umass_, const range &rgeps_,
      const range &rgrho_, const range &rgye_);

  // destructor
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline ~eos_idealgas();

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_valid_rho_eps_ye(
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL eps, ///< Specific internal energy \f$ \epsilon \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_press_ye(
      const CCTK_REAL rho,   ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL press, ///< Pressure \f$ P \f$
      const CCTK_REAL ye     ///< Electron fraction \f$ Y_e \f$
  ) const;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_valid_rho_eps_ye(
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL eps, ///< Specific internal energy \f$ \epsilon \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  temp_from_valid_rho_eps_ye(
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL eps, ///< Specific internal energy \f$ \epsilon \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  press_derivs_from_valid_rho_eps_ye(
      CCTK_REAL &press,  ///< Pressure \f$ P \f$
      CCTK_REAL &dpdrho, ///< Partial derivative \f$ \frac{\partial P}{\partial
                         ///< \rho} \f$
      CCTK_REAL &dpdeps, ///< Partial derivative \f$ \frac{\partial P}{\partial
                         ///< \epsilon} \f$
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL eps, ///< Specific internal energy \f$ \epsilon \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  entropy_from_valid_rho_temp_ye(
      const CCTK_REAL rho,  ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL temp, ///< Temperature \f$ T \f$
      const CCTK_REAL ye    ///< Electron fraction \f$ Y_e \f$
  ) const;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  entropy_from_valid_rho_eps_ye(
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL eps, ///< Specific internal energy \f$ \epsilon \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_temp_ye(
      const CCTK_REAL rho,  ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL temp, ///< Temperature \f$ T \f$ in MeV
      const CCTK_REAL ye    ///< Electron fraction \f$ Y_e \f$
  ) const;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline range
  range_eps_from_valid_rho_ye(
      const CCTK_REAL rho, ///< Rest mass density  \f$ \rho \f$
      const CCTK_REAL ye   ///< Electron fraction \f$ Y_e \f$
  ) const;
};

// constructor
CCTK_HOST
CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_idealgas::eos_idealgas(
    CCTK_REAL gamma_, CCTK_REAL umass_, const range &rgeps_,
    const range &rgrho_, const range &rgye_)
    : gamma(gamma_), gm1(gamma_ - 1), rgeps(rgeps_) {
  if (gamma < 1) {
    assert(0);
    // runtime_error("EOS_IdealGas: initialized with gamma < 1");
  }
  if (gamma > 2) { // Ensure subluminal Soundspeed and P < E
    rgeps.max = min(rgeps.max, 1 / (gamma * (gamma - 2)));
  }
  // set_range_h(range(1 + gamma * rgeps.min, 1 + gamma * rgeps.max));
  set_range_rho(rgrho_);
  set_range_ye(rgye_);
  temp_over_eps = gm1 * umass_;
  set_range_temp(range(temp_over_eps * rgeps.min, temp_over_eps * rgeps.max));
}

// destructor
CCTK_HOST CCTK_DEVICE
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_idealgas::~eos_idealgas() {}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_idealgas::press_from_valid_rho_eps_ye(const CCTK_REAL rho,
                                          const CCTK_REAL eps,
                                          const CCTK_REAL ye) const {
  return gm1 * rho * eps;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_idealgas::eps_from_valid_rho_press_ye(const CCTK_REAL rho,
                                          const CCTK_REAL press,
                                          const CCTK_REAL ye) const {
  return press / (rho * gm1);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_idealgas::csnd_from_valid_rho_eps_ye(const CCTK_REAL rho,
                                         const CCTK_REAL eps,
                                         const CCTK_REAL ye) const {
  return sqrt(gm1 * eps / (eps + 1 / gamma));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_idealgas::temp_from_valid_rho_eps_ye(const CCTK_REAL rho,
                                         const CCTK_REAL eps,
                                         const CCTK_REAL ye) const {
  return temp_over_eps * eps;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
eos_idealgas::press_derivs_from_valid_rho_eps_ye(
    CCTK_REAL &press, CCTK_REAL &dpdrho, CCTK_REAL &dpdeps, const CCTK_REAL rho,
    const CCTK_REAL eps, const CCTK_REAL ye) const {
  press = gm1 * rho * eps;
  dpdrho = gm1 * eps;
  dpdeps = gm1 * rho;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_idealgas::entropy_from_valid_rho_temp_ye(const CCTK_REAL rho,
                                             const CCTK_REAL temp,
                                             const CCTK_REAL ye) const {
  return log(temp * pow(rho, -gm1) / temp_over_eps);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_idealgas::entropy_from_valid_rho_eps_ye(const CCTK_REAL rho,
                                            const CCTK_REAL eps,
                                            const CCTK_REAL ye) const {
  return log(eps * pow(rho, -gm1));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_idealgas::eps_from_valid_rho_temp_ye(const CCTK_REAL rho,
                                         const CCTK_REAL temp,
                                         const CCTK_REAL ye) const {
  return temp / temp_over_eps;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos::range
eos_idealgas::range_eps_from_valid_rho_ye(const CCTK_REAL rho,
                                          const CCTK_REAL ye) const {
  return rgeps;
}

} // namespace EOSX

#endif
