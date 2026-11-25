#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

#include "atmo_global.hxx"
#include "setup_eos.hxx"

namespace AsterX {
using namespace EOSX;
using namespace Con2PrimFactory;

atmo_global global_atmo;

namespace {

void CheckAtmoParams() {
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(evolution_eos, "Tabulated3d") &&
      (!thermal_eos_atmo || use_press_atmo))
    CCTK_ERROR("Tabulated3d atmosphere requires thermal_eos_atmo=yes "
               "and use_press_atmo=no");
  if (Ye_atmo_beq && !CCTK_EQUALS(evolution_eos, "Tabulated3d"))
    CCTK_ERROR("Ye_atmo_beq is supported only for Tabulated3d");

  if (!std::isfinite(rho_abs_min) || rho_abs_min < 0.0 ||
      (!Ye_atmo_beq && !std::isfinite(Ye_atmo)) ||
      !std::isfinite(atmo_tol) || atmo_tol < 0.0 ||
      !std::isfinite(r_atmo) || r_atmo <= 0.0 ||
      !std::isfinite(n_rho_atmo) || n_rho_atmo < 0.0)
    CCTK_ERROR("Invalid atmosphere density, composition, cutoff or grading");

  if (thermal_eos_atmo) {
    const CCTK_REAL value = use_press_atmo ? p_atmo : t_atmo;
    const CCTK_REAL exponent = use_press_atmo ? n_press_atmo : n_temp_atmo;
    if (!std::isfinite(value) || value < 0.0 ||
        !std::isfinite(exponent) || exponent < 0.0)
      CCTK_ERROR("Invalid active atmosphere pressure/temperature or grading");
  } else if (CCTK_EQUALS(initial_data_eos, "Polytropic")) {
    if (!std::isfinite(poly_gamma) || poly_gamma <= 1.0 ||
        !std::isfinite(poly_k) || poly_k <= 0.0)
      CCTK_ERROR("Cold polytropic atmosphere requires gamma > 1 and K > 0");
  }

  if (global_atmo.valid &&
      !std::isfinite(global_atmo.atmo.rho_atmo * (1 + atmo_tol)))
    CCTK_ERROR("Atmosphere density cutoff is not finite");
}

template <typename EOSIDType, typename EOSType>
void SetupAtmo(const EOSIDType *eos_1p, const EOSType *eos_3p) {
  DECLARE_CCTK_PARAMETERS;

  if (!eos_3p || (!thermal_eos_atmo && !eos_1p))
    CCTK_ERROR("Required EOS was not initialized for atmosphere construction");

  const auto valid_range = [](const auto &rg) {
    return std::isfinite(rg.min) && std::isfinite(rg.max) && rg.min <= rg.max;
  };
  if (!valid_range(eos_3p->rgrho) || eos_3p->rgrho.min <= 0.0 ||
      !valid_range(eos_3p->rgtemp) || eos_3p->rgtemp.min < 0.0 ||
      !valid_range(eos_3p->rgye) || !valid_range(eos_3p->rgeps))
    CCTK_ERROR("Invalid EOS bounds for atmosphere construction");

  const auto atmo = make_atmo(
      eos_1p, eos_3p, 0.0, rho_abs_min, p_atmo, t_atmo, Ye_atmo, r_atmo,
      n_rho_atmo, n_press_atmo, n_temp_atmo, atmo_tol, thermal_eos_atmo,
      use_press_atmo, Ye_atmo_beq);
  if (!std::isfinite(atmo.rho_atmo) || atmo.rho_atmo <= 0.0 ||
      !std::isfinite(atmo.eps_atmo) ||
      !std::isfinite(atmo.press_atmo) || atmo.press_atmo < 0.0 ||
      !std::isfinite(atmo.temp_atmo) || atmo.temp_atmo < 0.0 ||
      !std::isfinite(atmo.ye_atmo) || !std::isfinite(atmo.entropy_atmo) ||
      !std::isfinite(atmo.rho_cut))
    CCTK_ERROR("Atmosphere construction did not produce a finite EOS state");

  global_atmo.store(eos_1p, eos_3p, atmo, rho_abs_min, p_atmo, t_atmo, Ye_atmo,
                     thermal_eos_atmo, use_press_atmo, Ye_atmo_beq);
  CCTK_VINFO("Prepared inner atmosphere: rho=%.16e eps=%.16e P=%.16e "
             "T=%.16e Ye=%.16e kappa=%.16e",
             atmo.rho_atmo, atmo.eps_atmo, atmo.press_atmo,
             atmo.temp_atmo, atmo.ye_atmo, atmo.entropy_atmo);
}

} // namespace

extern "C" void AsterX_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  // CarpetX runs PARAMCHECK after EOS setup on fresh starts and recovery.
  global_atmo.valid = false;
  if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    if (Ye_atmo_beq)
      CCTK_ERROR("Ye_atmo_beq is supported only for Tabulated3d");
    return; // Hybrid temperature operations are not implemented.
  }

  CheckAtmoParams();
  const auto setup = [](const auto *eos_3p) {
    if (global_eos_1p_pwpoly)
      SetupAtmo(global_eos_1p_pwpoly, eos_3p);
    else
      SetupAtmo(global_eos_1p_poly, eos_3p);
  };
  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    if (!global_eos_3p_ig ||
        !std::isfinite(global_eos_3p_ig->temp_over_eps) ||
        global_eos_3p_ig->temp_over_eps <= 0.0)
      CCTK_ERROR("Ideal-gas atmosphere requires a positive finite T/eps");
    setup(global_eos_3p_ig);
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    if (!global_eos_3p_tab3d || !(global_eos_3p_tab3d->rgtemp.min > 0.0))
      CCTK_ERROR("Tabulated atmosphere requires a positive minimum temperature");
    setup(global_eos_3p_tab3d);
  } else {
    CCTK_ERROR("Unsupported evolution EOS for atmosphere construction");
  }
}

bool get_global_atmo(const void *eos_1p, const void *eos_3p,
                     atmosphere &atmo) {
  DECLARE_CCTK_PARAMETERS;

  if (global_atmo.load(eos_1p, eos_3p, rho_abs_min, p_atmo, t_atmo, Ye_atmo,
                        n_rho_atmo, n_press_atmo, n_temp_atmo, atmo_tol,
                        thermal_eos_atmo, use_press_atmo, Ye_atmo_beq, atmo))
    return true;

  // Recheck changed/graded settings before the caller uses the device builder.
  // No EOS table is accessed here, and the shared state is never updated.
  if (!CCTK_EQUALS(evolution_eos, "Hybrid"))
    CheckAtmoParams();
  return false;
}

} // namespace AsterX
