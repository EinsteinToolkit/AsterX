#ifndef ASTERX_ATMO_GLOBAL_HXX
#define ASTERX_ATMO_GLOBAL_HXX

#include "atmo.hxx"

namespace AsterX {

// Filled during PARAMCHECK; grid calls only read and copy this state.
struct atmo_global {
  Con2PrimFactory::atmosphere atmo{};
  const void *cold = nullptr;
  const void *eos = nullptr;
  CCTK_REAL rho = 0.0;
  CCTK_REAL press = 0.0;
  CCTK_REAL temp = 0.0;
  CCTK_REAL Ye = 0.0;
  bool thermal = false;
  bool use_press = false;
  bool Ye_beq = false;
  bool valid = false;

  void store(const void *eos_1p, const void *eos_3p,
             const Con2PrimFactory::atmosphere &state,
             CCTK_REAL rho_abs_min, CCTK_REAL p_atmo,
             CCTK_REAL t_atmo, CCTK_REAL Ye_atmo,
             bool thermal_eos_atmo, bool use_press_atmo,
             bool Ye_atmo_beq) {
    atmo = state;
    cold = eos_1p;
    eos = eos_3p;
    rho = rho_abs_min;
    press = p_atmo;
    temp = t_atmo;
    Ye = Ye_atmo;
    thermal = thermal_eos_atmo;
    use_press = use_press_atmo;
    Ye_beq = Ye_atmo_beq;
    valid = true;
  }

  bool load(const void *eos_1p, const void *eos_3p,
            CCTK_REAL rho_abs_min, CCTK_REAL p_atmo,
            CCTK_REAL t_atmo, CCTK_REAL Ye_atmo,
            CCTK_REAL n_rho_atmo, CCTK_REAL n_press_atmo,
            CCTK_REAL n_temp_atmo, CCTK_REAL atmo_tol,
            bool thermal_eos_atmo, bool use_press_atmo, bool Ye_atmo_beq,
            Con2PrimFactory::atmosphere &state) const {
    if (!valid || eos != eos_3p || rho != rho_abs_min ||
        Ye_beq != Ye_atmo_beq || (!Ye_beq && Ye != Ye_atmo) ||
        thermal != thermal_eos_atmo || n_rho_atmo != 0.0)
      return false;

    // Only active inputs and exponents affect the atmosphere.
    if (thermal) {
      if (use_press != use_press_atmo ||
          (use_press ? press != p_atmo || n_press_atmo != 0.0
                     : temp != t_atmo || n_temp_atmo != 0.0))
        return false;
    } else if (cold != eos_1p) {
      return false;
    }

    const CCTK_REAL rho_cut = atmo.rho_atmo * (1 + atmo_tol);
    if (!std::isfinite(atmo_tol) || atmo_tol < 0.0 ||
        !std::isfinite(rho_cut))
      return false;
    state = atmo;
    // Steering the tolerance changes the copy, not the stored thermodynamics.
    state.rho_cut = rho_cut;
    return true;
  }
};

extern atmo_global global_atmo;

bool get_global_atmo(const void *eos_1p, const void *eos_3p,
                     Con2PrimFactory::atmosphere &atmo);

} // namespace AsterX

#endif
