#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>

#include <setup_eos.hxx>

namespace AsterSeeds {
using namespace Loop;
using namespace EOSX;

enum class eos_3param { IdealGas, Hybrid, Tabulated };

template <typename EOSType>
void AsterSeeds_SetTemp_typeEoS(CCTK_ARGUMENTS, EOSType *eos_3p) {

  DECLARE_CCTK_ARGUMENTSX_SetTemp;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL epsL = eps(p.I);
        temperature(p.I) =
            eos_3p->temp_from_rho_eps_ye(rho(p.I), epsL, Ye(p.I));
        eps(p.I) = epsL;
      });
}

template <typename EOSType>
void AsterSeeds_SetEntropy_typeEoS(CCTK_ARGUMENTS, EOSType *eos_3p) {

  DECLARE_CCTK_ARGUMENTSX_SetEntropy;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // kappa ~ eps * rho^(1-gamma) is 0*inf = NaN in vacuum (e.g. outside
        // a TOV star); those points are reset by the atmosphere later
        if (rho(p.I) > 0) {
          CCTK_REAL epsL = eps(p.I);
          entropy(p.I) =
              eos_3p->kappa_from_rho_eps_ye(rho(p.I), epsL, Ye(p.I));
          eps(p.I) = epsL;
        } else {
          entropy(p.I) = 0;
        }
      });
}

extern "C" void SetYe(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTSX_SetYe;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { Ye(p.I) = Ye_atmo; });
}

extern "C" void SetTemp(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTSX_SetTemp;
  DECLARE_CCTK_PARAMETERS;

  // defining EOS objects
  eos_3param eos_3p_type;

  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eos_3p_type = eos_3param::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eos_3p_type = eos_3param::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    eos_3p_type = eos_3param::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  switch (eos_3p_type) {
  case eos_3param::IdealGas: {
    auto eos_3p_ig = global_eos_3p_ig;
    AsterSeeds_SetTemp_typeEoS(CCTK_PASS_CTOC, eos_3p_ig);
    break;
  }
  case eos_3param::Hybrid: {
    if (global_eos_3p_hyb_pwpoly) {
      auto eos_3p_hyb = global_eos_3p_hyb_pwpoly;
      AsterSeeds_SetTemp_typeEoS(CCTK_PASS_CTOC, eos_3p_hyb);
    } else if (global_eos_3p_hyb_poly) {
      auto eos_3p_hyb = global_eos_3p_hyb_poly;
      AsterSeeds_SetTemp_typeEoS(CCTK_PASS_CTOC, eos_3p_hyb);
    } else {
      CCTK_ERROR(
          "Hybrid EOS selected but no hybrid EOS object was initialized");
    }
    break;
  }
  case eos_3param::Tabulated: {
    auto eos_3p_tab3d = global_eos_3p_tab3d;
    AsterSeeds_SetTemp_typeEoS(CCTK_PASS_CTOC, eos_3p_tab3d);
    break;
  }
  default:
    assert(0);
  }
}

extern "C" void SetEntropy(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTSX_SetEntropy;
  DECLARE_CCTK_PARAMETERS;

  // defining EOS objects
  eos_3param eos_3p_type;

  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eos_3p_type = eos_3param::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eos_3p_type = eos_3param::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    eos_3p_type = eos_3param::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  switch (eos_3p_type) {
  case eos_3param::IdealGas: {
    auto eos_3p_ig = global_eos_3p_ig;
    AsterSeeds_SetEntropy_typeEoS(CCTK_PASS_CTOC, eos_3p_ig);
    break;
  }
  case eos_3param::Hybrid: {
    if (global_eos_3p_hyb_pwpoly) {
      auto eos_3p_hyb = global_eos_3p_hyb_pwpoly;
      AsterSeeds_SetEntropy_typeEoS(CCTK_PASS_CTOC, eos_3p_hyb);
    } else if (global_eos_3p_hyb_poly) {
      auto eos_3p_hyb = global_eos_3p_hyb_poly;
      AsterSeeds_SetEntropy_typeEoS(CCTK_PASS_CTOC, eos_3p_hyb);
    } else {
      CCTK_ERROR(
          "Hybrid EOS selected but no hybrid EOS object was initialized");
    }
    break;
  }
  case eos_3param::Tabulated: {
    auto eos_3p_tab3d = global_eos_3p_tab3d;
    AsterSeeds_SetEntropy_typeEoS(CCTK_PASS_CTOC, eos_3p_tab3d);
    break;
  }
  default:
    assert(0);
  }
}

} // namespace AsterSeeds
