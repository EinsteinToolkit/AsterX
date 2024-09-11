#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

#include "seeds_utils.hxx"

namespace AsterSeeds {
using namespace std;
using namespace Loop;
using namespace AsterUtils;

extern "C" void AsterSeeds_InitializeCenteredAvec_TOV(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_InitializeCenteredAvec_TOV;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(Afield_config, "internal dipole")) {

    /* computing cell centered vector potential components */
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL x_local = p.x - dipole_x[0];
          CCTK_REAL y_local = p.y - dipole_y[0];
          CCTK_REAL Pcut = press_max * press_cut;
          CCTK_REAL Pdiff = std::max(press(p.I) - Pcut, 0.0);
          CCTK_REAL Aphi_local = Ab * pow(Pdiff, Avec_kappa);
          Avec_x_cent(p.I) = -y_local * Aphi_local;
          Avec_y_cent(p.I) = x_local * Aphi_local;
          Avec_z_cent(p.I) = 0.0;
        });

  } else if (CCTK_EQUALS(Afield_config, "external dipole")) {

    /* computing cell centered vector potential components */
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL x_local = p.x - dipole_x[0];
          CCTK_REAL y_local = p.y - dipole_y[0];
          CCTK_REAL z_local = p.z - dipole_z[0];
          CCTK_REAL cylrad2 = x_local * x_local + y_local * y_local;
          CCTK_REAL rsph =
              sqrt(x_local * x_local + y_local * y_local + z_local * z_local);
          CCTK_REAL rsph3 = pow(rsph, 3.0);
          CCTK_REAL r03 = pow(r0, 3.0);
          CCTK_REAL Aphi_local =
              B0 * (r03 / (r03 + rsph3)) / sqrt(cylrad2 + 1.0e-16);
          Avec_x_cent(p.I) = -y_local * Aphi_local;
          Avec_y_cent(p.I) = x_local * Aphi_local;
          Avec_z_cent(p.I) = 0.0;
        });

  } else {
    CCTK_ERROR("Vector potential configuration not defined");
  }
}

extern "C" void AsterSeeds_InitializeStagAvec_TOV(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_InitializeStagAvec_TOV;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_x(p.I) = calc_avg_c2e(Avec_x_cent, p, 0);
      });

  grid.loop_int<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_y(p.I) = calc_avg_c2e(Avec_y_cent, p, 1);
      });

  grid.loop_int<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_z(p.I) = calc_avg_c2e(Avec_z_cent, p, 2);
      });
}

} // namespace AsterSeeds
