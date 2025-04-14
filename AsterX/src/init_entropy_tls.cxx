#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

#include "aster_utils.hxx"

namespace AsterX {
using namespace Loop;
using namespace std;

enum class eos_t { IdealGas, Hybrid, Tabulated };

extern "C" void AsterX_InitEntropyTLs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_InitEntropyTLs;
  DECLARE_CCTK_PARAMETERS;

  eos_t eostype;
  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eostype = eos_t::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eostype = eos_t::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated")) {
    eostype = eos_t::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        switch (eostype) {
        case eos_t::IdealGas: {
          phys_ent(p.I) = std::log(entropy(p.I));
          break;
        }
        case eos_t::Hybrid: {
          printf("Hybrid EOS is not yet supported");
          assert(0);
          break;
        }
        case eos_t::Tabulated: {
          printf("Tabulated EOS is not yet supported");
          assert(0);
          break;
        }
        default:
          assert(0);
        }

        // Simple copy
        ent_m1(p.I) = phys_ent(p.I);
        ent_m2(p.I) = phys_ent(p.I);

        // Init Residual to 0
        r_ent(p.I) = 0.0;
      });
}

} // namespace AsterX
