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

extern "C" void AsterSeeds_InitializeCenteredAvec_BNS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_InitializeCenteredAvec_BNS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(Afield_config, "internal dipole")) {

    /* computing cell centered vector potential components */
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          // dummy initialization
          CCTK_REAL x_local = 0.0;
          CCTK_REAL y_local = 0.0;

          // For star 1 at minus side
          if (p.x < 0) {
            x_local = p.x - dipole_x[0];
            y_local = p.y - dipole_y[0];
          }
          // For star 2 at minus side
          else {
            x_local = p.x - dipole_x[1];
            y_local = p.y - dipole_y[1];
          }

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
          // For star 1 at minus side
          CCTK_REAL x_local_s1 = p.x - dipole_x[0];
          CCTK_REAL y_local_s1 = p.y - dipole_y[0];
          CCTK_REAL z_local_s1 = p.z - dipole_z[0];
          CCTK_REAL cylrad2_s1 =
              x_local_s1 * x_local_s1 + y_local_s1 * y_local_s1;
          CCTK_REAL rsph_s1 =
              sqrt(x_local_s1 * x_local_s1 + y_local_s1 * y_local_s1 +
                   z_local_s1 * z_local_s1);

          CCTK_REAL Aphi_local_s1 =
              B0 * (pow(r0, 3.0) / (pow(r0, 3.0) + pow(rsph_s1, 3.0))) /
              sqrt(cylrad2_s1 + 1.0e-16);

          // For star 2 at minus side
          CCTK_REAL x_local_s2 = p.x - dipole_x[1];
          CCTK_REAL y_local_s2 = p.y - dipole_y[1];
          CCTK_REAL z_local_s2 = p.z - dipole_z[1];
          CCTK_REAL cylrad2_s2 =
              x_local_s2 * x_local_s2 + y_local_s2 * y_local_s2;
          CCTK_REAL rsph_s2 =
              sqrt(x_local_s2 * x_local_s2 + y_local_s2 * y_local_s2 +
                   z_local_s2 * z_local_s2);

          CCTK_REAL Aphi_local_s2 =
              B0 * (pow(r0, 3.0) / (pow(r0, 3.0) + pow(rsph_s2, 3.0))) /
              sqrt(cylrad2_s2 + 1.0e-16);

          Avec_x_cent(p.I) =
              -(y_local_s1 * Aphi_local_s1 + y_local_s2 * Aphi_local_s2);
          Avec_y_cent(p.I) =
              x_local_s1 * Aphi_local_s1 + x_local_s2 * Aphi_local_s2;
          Avec_z_cent(p.I) = 0.0;
        });

  } else {
    CCTK_ERROR("Vector potential configuration not defined");
  }
}

extern "C" void AsterSeeds_InitializeStagAvec_BNS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_InitializeStagAvec_BNS;
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
