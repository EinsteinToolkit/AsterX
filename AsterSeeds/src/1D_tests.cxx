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

extern "C" void Tests1D_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Tests1D_Initialize;
  DECLARE_CCTK_PARAMETERS;

  // For all the tests, the initial data EOS is ideal gas
  // Constructing the IG EOS object
  eos::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
      rgye(ye_min, ye_max);

  const eos_idealgas eos_th(gl_gamma, particle_mass, rgeps, rgrho, rgye);
  const CCTK_REAL dummy_ye = 0.5;

  if (CCTK_EQUALS(test_case, "equilibrium")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0;
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
          press(p.I) = 1.0;
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

  } else if (CCTK_EQUALS(test_case, "sound wave")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0 + amplitude * sin(M_PI * p.x);
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
          press(p.I) = 1.0; // should add kinetic energy here
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

  } else if (CCTK_EQUALS(test_case, "Alfven wave")) {
    const CCTK_REAL A0 = 1.0;
    const CCTK_REAL va = 0.5;
    const CCTK_REAL k = 2 * M_PI;

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0;
          vely(p.I) = -va * A0 * cos(k * p.x);
          velz(p.I) = -va * A0 * sin(k * p.x);
          press(p.I) = 0.5; // should add kinetic energy here
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          Avec_x(p.I) = p.z * cos(k * p.x) - p.y * sin(k * p.x);
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = -p.z / 2.0; });
    // CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = p.x*sin(k*p.x) - p.z; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = p.y / 2.0; });
    // CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = p.y - p.x*cos(k*p.x); });

  } else if (CCTK_EQUALS(test_case, "shock tube")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = 2.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 2.0;
          } else {
            rho(p.I) = 1.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 1.0;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

  } else if (CCTK_EQUALS(test_case, "Balsara1") ||
             CCTK_EQUALS(test_case, "Balsara2") ||
             CCTK_EQUALS(test_case, "Balsara3") ||
             CCTK_EQUALS(test_case, "Balsara4") ||
             CCTK_EQUALS(test_case, "Balsara5")) {

    auto myNaN = numeric_limits<CCTK_REAL>::quiet_NaN();

    CCTK_REAL rhol = myNaN;
    CCTK_REAL vxl = myNaN;
    CCTK_REAL vyl = myNaN;
    CCTK_REAL vzl = myNaN;
    CCTK_REAL pressl = myNaN;
    CCTK_REAL Bxl = myNaN;
    CCTK_REAL Byl = myNaN;
    CCTK_REAL Bzl = myNaN;

    CCTK_REAL rhor = myNaN;
    CCTK_REAL vxr = myNaN;
    CCTK_REAL vyr = myNaN;
    CCTK_REAL vzr = myNaN;
    CCTK_REAL pressr = myNaN;
    CCTK_REAL Bxr = myNaN;
    CCTK_REAL Byr = myNaN;
    CCTK_REAL Bzr = myNaN;

    if (CCTK_EQUALS(test_case, "Balsara1")) {
      rhol = 1.0;
      vxl = 0.0;
      vyl = 0.0;
      vzl = 0.0;
      pressl = 1.0;
      Bxl = 0.5;
      Byl = 1.0;
      Bzl = 0.0;

      rhor = 0.125;
      vxr = 0.0;
      vyr = 0.0;
      vzr = 0.0;
      pressr = 0.1;
      Bxr = 0.5;
      Byr = -1.0;
      Bzr = 0.0;
    } else if (CCTK_EQUALS(test_case, "Balsara2")) {
      rhol = 1.0;
      vxl = 0.0;
      vyl = 0.0;
      vzl = 0.0;
      pressl = 30.0;
      Bxl = 5.0;
      Byl = 6.0;
      Bzl = 6.0;

      rhor = 1.0;
      vxr = 0.0;
      vyr = 0.0;
      vzr = 0.0;
      pressr = 1.0;
      Bxr = 5.0;
      Byr = 0.7;
      Bzr = 0.7;
    } else if (CCTK_EQUALS(test_case, "Balsara3")) {
      rhol = 1.0;
      vxl = 0.0;
      vyl = 0.0;
      vzl = 0.0;
      pressl = 1000.0;
      Bxl = 10.0;
      Byl = 7.0;
      Bzl = 7.0;

      rhor = 1.0;
      vxr = 0.0;
      vyr = 0.0;
      vzr = 0.0;
      pressr = 0.1;
      Bxr = 10.0;
      Byr = 0.7;
      Bzr = 0.7;
    } else if (CCTK_EQUALS(test_case, "Balsara4")) {
      rhol = 1.0;
      vxl = 0.999;
      vyl = 0.0;
      vzl = 0.0;
      pressl = 0.1;
      Bxl = 10.0;
      Byl = 7.0;
      Bzl = 7.0;

      rhor = 1.0;
      vxr = -0.999;
      vyr = 0.0;
      vzr = 0.0;
      pressr = 0.1;
      Bxr = 10.0;
      Byr = -7.0;
      Bzr = -7.0;
    } else if (CCTK_EQUALS(test_case, "Balsara5")) {
      rhol = 1.08;
      vxl = 0.4;
      vyl = 0.3;
      vzl = 0.2;
      pressl = 0.95;
      Bxl = 2.0;
      Byl = 0.3;
      Bzl = 0.3;

      rhor = 1.0;
      vxr = -0.45;
      vyr = -0.2;
      vzr = 0.2;
      pressr = 1.0;
      Bxr = 2.0;
      Byr = -0.7;
      Bzr = 0.5;
    } else {
      CCTK_ERROR("Balsara type not defined");
    }

    CCTK_REAL theta_x = myNaN;
    CCTK_REAL theta_y = myNaN;
    CCTK_REAL theta_z = myNaN;
    vector<int> rotate_seq{2, 1, 0};

    if (CCTK_EQUALS(shock_dir, "x")) {
      theta_z = 0.0;
      theta_y = 0.0;
      theta_x = 0.0;
    } else if (CCTK_EQUALS(shock_dir, "y")) {
      theta_z = 0.5 * M_PI;
      theta_y = 0.0;
      theta_x = 0.5 * M_PI;
    } else if (CCTK_EQUALS(shock_dir, "z")) {
      theta_z = -0.5 * M_PI;
      theta_y = -0.5 * M_PI;
      theta_x = 0.0;
    } else {
      theta_z = rotate_angle_z * M_PI;
      theta_y = rotate_angle_y * M_PI;
      theta_x = rotate_angle_x * M_PI;
    }

    // Rotation matrix R_{ij} that rotates the coordinate system through a
    // counterclockwise angle of x0-, x1-, x2-axes. (v'_i = R_{ij}v_j under this
    // coordinate transformation).
    const auto calc_R = [=] CCTK_HOST CCTK_DEVICE(
                            const CCTK_REAL th_x, const CCTK_REAL th_y,
                            const CCTK_REAL th_z,
                            vector<int> seq) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      const mat<CCTK_REAL, 3> R_x{
          1.0, 0.0, 0.0, 0.0, cos(th_x), sin(th_x), 0.0, -sin(th_x), cos(th_x)};
      const mat<CCTK_REAL, 3> R_y{cos(th_y), 0.0, -sin(th_y), 0.0, 1.0, 0.0,
                                  sin(th_y), 0.0, cos(th_y)};
      const mat<CCTK_REAL, 3> R_z{
          cos(th_z), sin(th_z), 0.0, -sin(th_z), cos(th_z), 0.0, 0.0, 0.0, 1.0};
      const vec<mat<CCTK_REAL, 3>, 3> Rs{R_x, R_y, R_z};
      // R_{ij} = R3_{ik} R2_{kl} R1_{lj}
      return mat<CCTK_REAL, 3>([&](int i, int j) ARITH_INLINE {
        return sum<3>([&](int k) ARITH_INLINE {
          return Rs(seq[2])(i, k) * sum<3>([&](int l) ARITH_INLINE {
                   return Rs(seq[1])(k, l) * Rs(seq[0])(l, j);
                 });
        });
      });
    };

    // Tranform vector
    const auto transform_vec =
        [=] CCTK_HOST CCTK_DEVICE(const mat<CCTK_REAL, 3> R_mat,
                                  const vec<CCTK_REAL, 3> u_vec)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
              return vec<CCTK_REAL, 3>([&](int i) ARITH_INLINE {
                return sum<3>(
                    [&](int k) ARITH_INLINE { return R_mat(i, k) * u_vec(k); });
              });
            };

    // Construct rotation matrix and its inverse
    const auto R = calc_R(theta_x, theta_y, theta_z, rotate_seq);
    const auto Rinv =
        calc_R(-theta_x, -theta_y, -theta_z,
               vector<int>(rotate_seq.rbegin(), rotate_seq.rend()));

    // B and v in new coordinates: B' = Bs, v' = vs
    const vec<CCTK_REAL, 3> newBls{Bxl, Byl, Bzl};
    const vec<CCTK_REAL, 3> newBrs{Bxr, Byr, Bzr};
    const vec<CCTK_REAL, 3> newvls{vxl, vyl, vzl};
    const vec<CCTK_REAL, 3> newvrs{vxr, vyr, vzr};

    // B and v in old coordinates: B_i = Rinv_{ij}B'_j, v_i = Rinv_{ij}v'_j
    const auto oldBls = transform_vec(Rinv, newBls);
    const auto oldBrs = transform_vec(Rinv, newBrs);
    const auto oldvls = transform_vec(Rinv, newvls);
    const auto oldvrs = transform_vec(Rinv, newvrs);

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto newxs = transform_vec(R, vec<CCTK_REAL, 3>{p.x, p.y, p.z});

          if (newxs(0) <= 0.0) {
            rho(p.I) = rhol;
            velx(p.I) = oldvls(0);
            vely(p.I) = oldvls(1);
            velz(p.I) = oldvls(2);
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velx(p.I) = oldvrs(0);
            vely(p.I) = oldvrs(1);
            velz(p.I) = oldvrs(2);
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    // Set up Avec (only support oldBs that have at least one zero component):
    // we are set A in such a way that it only depends on the direction where
    // B_i=0 (no step). Then, when take derivatives of A, it's well-behaved.
    const CCTK_REAL tiny = 1e-12;
    const auto are_all_components_nonzero = [&](const vec<CCTK_REAL, 3> Bs) {
      return (abs(Bs(0)) > tiny) && (abs(Bs(1)) > tiny) && (abs(Bs(2)) > tiny);
    };

    if (are_all_components_nonzero(oldBls) ||
        are_all_components_nonzero(oldBrs)) {
      CCTK_ERROR("All non-zero B components not supported yet.");
    }

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto newxs = transform_vec(R, vec<CCTK_REAL, 3>{p.x, p.y, p.z});
          if (newxs(0) <= 0.0) {
            Avec_x(p.I) = (abs(oldBls(0)) < tiny)
                              ? 0.0
                              : oldBls(1) * (p.z) - oldBls(2) * (p.y);
          } else {
            Avec_x(p.I) = (abs(oldBrs(0)) < tiny)
                              ? 0.0
                              : oldBrs(1) * (p.z) - oldBrs(2) * (p.y);
          }
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto newxs = transform_vec(R, vec<CCTK_REAL, 3>{p.x, p.y, p.z});
          if (newxs(0) <= 0.0) {
            Avec_y(p.I) = (abs(oldBls(1)) < tiny)
                              ? 0.0
                              : oldBls(2) * (p.x) - oldBls(0) * (p.z);
          } else {
            Avec_y(p.I) = (abs(oldBrs(1)) < tiny)
                              ? 0.0
                              : oldBrs(2) * (p.x) - oldBrs(0) * (p.z);
          }
        });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto newxs = transform_vec(R, vec<CCTK_REAL, 3>{p.x, p.y, p.z});
          if (newxs(0) <= 0.0) {
            Avec_z(p.I) = (abs(oldBls(2)) < tiny)
                              ? 0.0
                              : (oldBls(0) * (p.y) - oldBls(1) * (p.x));
          } else {
            Avec_z(p.I) = (abs(oldBrs(2)) < tiny)
                              ? 0.0
                              : (oldBrs(0) * (p.y) - oldBrs(1) * (p.x));
          }
        });

  } else {
    CCTK_ERROR("Test case not defined");
  }
}

} // namespace AsterSeeds
