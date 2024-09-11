#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

#include "eos.hxx"
#include "eos_idealgas.hxx"
#include "seeds_utils.hxx"

namespace AsterSeeds {
using namespace std;
using namespace Loop;
using namespace EOSX;
using namespace AsterUtils;

extern "C" void Tests2D_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Tests2D_Initialize;
  DECLARE_CCTK_PARAMETERS;

  // For all the tests, the initial data EOS is ideal gas
  // Constructing the IG EOS object
  eos::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
      rgye(ye_min, ye_max);

  const eos_idealgas eos_th(gl_gamma, particle_mass, rgeps, rgrho, rgye);
  const CCTK_REAL dummy_ye = 0.5;

  // See Cipolletta et al (2020) and Del Zanna, Bucciantini, Londrillo (2003)
  if (CCTK_EQUALS(test_case, "magnetic rotor")) {
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const CCTK_REAL r_cyl = sqrt(pow2(p.x) + pow2(p.y));

          if (r_cyl < 0.1) {
            const CCTK_REAL omega = 9.95;
            rho(p.I) = 10.;
            velx(p.I) = -omega * p.y;
            vely(p.I) = omega * p.x;
            velz(p.I) = 0.;
          }

          else {
            rho(p.I) = 1.;
            velx(p.I) = 0.;
            vely(p.I) = 0.;
            velz(p.I) = 0.;
          }

          press(p.I) = 1.;
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.; });

    grid.loop_all<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.; });

    grid.loop_all<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = p.y; });
  }

  // See Cipolletta et al (2020), and Beckwith, Stone (2011)
  else if (CCTK_EQUALS(test_case, "magnetic loop advection")) {
    CCTK_REAL axial_vel;

    if (CCTK_EQUALS(mag_loop_adv_type, "2D")) {
      if (CCTK_EQUALS(mag_loop_adv_axial_vel, "zero"))
        axial_vel = 0.;
      else if (CCTK_EQUALS(mag_loop_adv_axial_vel, "non-zero"))
        axial_vel = 1. / 24.;
      else {
        CCTK_VERROR("Invalid loop advection case");
      }
    }

    // TODO: implement this (see Cipolletta et al (2020) and code)
    else if (CCTK_EQUALS(mag_loop_adv_type, "3D")) {
      CCTK_VERROR("Sorry, 3D loop advection hasn't been implemented yet");
    }

    else {
      CCTK_VERROR("Invalid loop advection type");
    }

    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.;
          press(p.I) = 3.;
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
          velx(p.I) = 1. / 12.0;
          vely(p.I) = 1. / 24.;
          velz(p.I) = axial_vel;
        });

    grid.loop_all<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.; });

    grid.loop_all<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.; });

    grid.loop_all<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const CCTK_REAL r_cyl = sqrt(pow2(p.x) + pow2(p.y));
          const CCTK_REAL R_loop = 0.3;
          const CCTK_REAL A_loop = 0.001;

          if (r_cyl < R_loop)
            Avec_z(p.I) = A_loop * (R_loop - r_cyl);
          else
            Avec_z(p.I) = 0.;
        });
  }

  // See Cipolletta et al (2020)
  else if (CCTK_EQUALS(test_case, "cylindrical blast")) {
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const CCTK_REAL r_cyl = sqrt(pow2(p.x) + pow2(p.y));
          const CCTK_REAL r_in = 0.8;
          const CCTK_REAL r_out = 1.0;
          const CCTK_REAL rho_in = 1e-2;
          const CCTK_REAL rho_out = 1e-4;
          const CCTK_REAL p_in = 1.0;
          const CCTK_REAL p_out = 3.0e-5;

          if (r_cyl < r_in) {
            rho(p.I) = rho_in;
            press(p.I) = p_in;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
          }

          else if (r_cyl > r_out) {
            rho(p.I) = rho_out;
            press(p.I) = p_out;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
          } else {
            rho(p.I) = exp(((r_out - r_cyl) * log(rho_in) +
                            (r_cyl - r_in) * log(rho_out)) /
                           (r_out - r_in));
            press(p.I) = exp(
                ((r_out - r_cyl) * log(p_in) + (r_cyl - r_in) * log(p_out)) /
                (r_out - r_in));
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.1 * p.y; });
  }

  else {
    CCTK_ERROR("Test case not defined");
  }
}

} // namespace AsterSeeds
