#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include "seeds_utils.hxx"

namespace AsterSeeds {
using namespace std;
using namespace Loop;

extern "C" void AsterSeeds_InitializeCenteredAvec(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_InitializeCenteredAvec;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(Afield_config, "internal dipole")) {

    /* computing cell centered vector potential components */
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL Pcut = press_max * press_cut;
          CCTK_REAL Pdiff = std::max(press(p.I) - Pcut, 0.0);
          CCTK_REAL Aphi_local = Ab * pow(Pdiff, Avec_kappa);
          Avec_x_cent(p.I) = -p.y * Aphi_local;
          Avec_y_cent(p.I) = p.x * Aphi_local;
          Avec_z_cent(p.I) = 0.0;
        });

  } else if (CCTK_EQUALS(Afield_config, "external dipole")) {

    /* computing cell centered vector potential components */
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL cylrad2 = p.x * p.x + p.y * p.y;
          CCTK_REAL rsph = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
          CCTK_REAL rsph3 = pow(rsph, 3.0);
          CCTK_REAL r03 = pow(r0, 3.0);
          CCTK_REAL Aphi_local =
              B0 * (r03 / (r03 + rsph3)) / sqrt(cylrad2 + 1.0e-16);
          Avec_x_cent(p.I) = -p.y * Aphi_local;
          Avec_y_cent(p.I) = p.x * Aphi_local;
          Avec_z_cent(p.I) = 0.0;
        });

  } else {
    CCTK_ERROR("Vector potential configuration not defined");
  }
}

extern "C" void AsterSeeds_InitializeStagAvec(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_InitializeStagAvec;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_x(p.I) = calc_avg_c2e(Avec_x_cent, p, 0);
      });

  grid.loop_int_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_y(p.I) = calc_avg_c2e(Avec_y_cent, p, 1);
      });

  grid.loop_int_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_z(p.I) = calc_avg_c2e(Avec_z_cent, p, 2);
      });
}

} // namespace AsterSeeds
