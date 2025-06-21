#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include "setup_eos.hxx"
#include "aster_utils.hxx"

namespace AsterX {
using namespace AsterUtils;
using namespace Loop;
using namespace EOSX;
using namespace std;

enum class eos_3param { IdealGas, Hybrid, Tabulated };

template <typename EOSIDType, typename EOSType>
void CheckPrims(CCTK_ARGUMENTS, EOSIDType *eos_1p, EOSType *eos_3p) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CheckPrims;
  DECLARE_CCTK_PARAMETERS;

  // Loop over the entire grid (0 to n-1 cells in each direction)
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        // Interpolate metric terms from vertices to center
        const smat<CCTK_REAL, 3> g{calc_avg_v2c(gxx, p), calc_avg_v2c(gxy, p),
                                   calc_avg_v2c(gxz, p), calc_avg_v2c(gyy, p),
                                   calc_avg_v2c(gyz, p), calc_avg_v2c(gzz, p)};

        vec<CCTK_REAL, 3> v_up{velx(p.I), vely(p.I), velz(p.I)};

        CCTK_REAL rhoL = rho(p.I);
        CCTK_REAL epsL = eps(p.I);
        CCTK_REAL pressL = press(p.I);
        CCTK_REAL entropyL = entropy(p.I);
        CCTK_REAL YeL = Ye(p.I);
        CCTK_REAL tempL = temperature(p.I);

        // Setting up atmosphere
        CCTK_REAL rho_atm = 0.0;   // dummy initialization
        CCTK_REAL press_atm = 0.0; // dummy initialization
        CCTK_REAL eps_atm = 0.0;   // dummy initialization
        CCTK_REAL temp_atm = 0.0;  // dummy initialization
        auto rgeps = eos_3p->range_eps_from_valid_rho_ye(rhoL, YeL);

        CCTK_REAL radial_distance = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);

        // Grading rho
        rho_atm =
            (radial_distance > r_atmo)
                ? (rho_abs_min * pow((r_atmo / radial_distance), n_rho_atmo))
                : rho_abs_min;
        rho_atm = std::max(eos_3p->rgrho.min, rho_atm);

        // Grading pressure based on either cold or thermal EOS
        if (thermal_eos_atmo) {
          // rho_atm = max(rho_atm, eos_3p->interptable->xmin<0>());
          temp_atm = (radial_distance > r_atmo)
                         ? (t_atmo * pow(r_atmo / radial_distance, n_temp_atmo))
                         : t_atmo;
          temp_atm = std::max(eos_3p->rgtemp.min, temp_atm);
          // temp_atm = max(temp_atm, eos_3p->interptable->xmin<1>());
          press_atm =
              eos_3p->press_from_valid_rho_temp_ye(rho_atm, temp_atm, Ye_atmo);
          eps_atm =
              eos_3p->eps_from_valid_rho_temp_ye(rho_atm, temp_atm, Ye_atmo);
          // eps_atm should be kept consistent with temp_atm, so we do not use
          // the setting below
          // eps_atm =
          //    std::min(std::max(eos_3p->rgeps.min, eps_atm),
          //    eos_3p->rgeps.max);
        } else {
          const CCTK_REAL gm1 = eos_1p->gm1_from_valid_rho(rho_atm);
          eps_atm = eos_1p->sed_from_valid_gm1(gm1);
          eps_atm =
              std::min(std::max(eos_3p->rgeps.min, eps_atm), eos_3p->rgeps.max);
          press_atm =
              eos_3p->press_from_valid_rho_eps_ye(rho_atm, eps_atm, Ye_atmo);
        }
        const CCTK_REAL rho_atmo_cut = rho_atm * (1 + atmo_tol);

        CCTK_REAL rhomax = eos_3p->rgrho.max;
        CCTK_REAL tempmax = eos_3p->rgtemp.max;
        CCTK_REAL epsmax = rgeps.max;


        // Lower velocity
        vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);

        // Check validity of primitives ----------
        // Mostly copy of code from c2p::prims_floors_and_ceilings

        const CCTK_REAL w_lim = sqrt(1.0 + vw_lim * vw_lim);
        const CCTK_REAL v_lim = vw_lim / w_lim;

        // ----------
        // Floor and ceiling for rho and velocity
        // Keeps pressure the same and changes eps
        // ----------

        // check if computed velocities are within the specified limit
        CCTK_REAL vsq_Sol = calc_contraction(v_low, v_up);
        CCTK_REAL sol_v = sqrt(vsq_Sol);
        if (sol_v > v_lim) {

          v_up *= v_lim / sol_v;
          v_low = v_low * v_lim / sol_v;
        }

        if (rhoL > rhomax) {

          // remove mass
          rhoL = rhomax;

          if (use_temperature) {
            epsL = eos_3p->eps_from_valid_rho_temp_ye(rhoL, tempL, YeL);
            pressL = eos_3p->press_from_valid_rho_temp_ye(rhoL, tempL, YeL);
          } else {
            epsL = eos_3p->eps_from_valid_rho_press_ye(rhoL, pressL, YeL);
          }
          entropyL = eos_3p->kappa_from_valid_rho_eps_ye(rhoL, epsL, YeL);
        }

        if (rhoL < rho_atmo_cut) {

          // add mass
          rhoL = rho_atm;

          if (use_temperature) {
            epsL = eos_3p->eps_from_valid_rho_temp_ye(rhoL, tempL, YeL);
            pressL = eos_3p->press_from_valid_rho_temp_ye(rhoL, tempL, YeL);
          } else {
            epsL = eos_3p->eps_from_valid_rho_press_ye(rhoL, pressL, YeL);
          }
          entropyL = eos_3p->kappa_from_valid_rho_eps_ye(rhoL, epsL, YeL);
        }

        // ----------
        // Floor and ceiling for temp
        // Keeps rho the same and changes press
        // ----------

        if (use_temperature) {
          // check the validity of the computed temperature
          if (tempL > tempmax) {
            tempL = tempmax;
            epsL = eos_3p->eps_from_valid_rho_temp_ye(rhoL, tempL, YeL);
            pressL = eos_3p->press_from_valid_rho_temp_ye(rhoL, tempL, YeL);
            entropyL = eos_3p->kappa_from_valid_rho_eps_ye(rhoL, epsL, YeL);
          }
          if (tempL < temp_atm) {
            tempL = temp_atm;
            epsL = eos_3p->eps_from_valid_rho_temp_ye(rhoL, tempL, YeL);
            pressL = eos_3p->press_from_valid_rho_temp_ye(rhoL, tempL, YeL);
            entropyL = eos_3p->kappa_from_valid_rho_eps_ye(rhoL, epsL, YeL);
          }
        }

        // ----------
        // Floor and ceiling for eps
        // Keeps rho the same and changes press
        // ----------

        // check the validity of the computed eps
        if (epsL > epsmax) {

          epsL = epsmax;
          tempL = eos_3p->temp_from_valid_rho_eps_ye(rhoL, epsL, YeL);
          pressL = eos_3p->press_from_valid_rho_eps_ye(rhoL, epsL, YeL);
          entropyL = eos_3p->kappa_from_valid_rho_eps_ye(rhoL, epsL, YeL);

        } else if (epsL < eps_atm) {

          epsL = eps_atm;
          tempL = eos_3p->temp_from_valid_rho_eps_ye(rhoL, epsL, YeL);
          pressL = eos_3p->press_from_valid_rho_eps_ye(rhoL, epsL, YeL);
          entropyL = eos_3p->kappa_from_valid_rho_eps_ye(rhoL, epsL, YeL);
        }

        // ----------
        // Floor and ceiling for Ye
        // ----------

        // TODO:

        // ---------- End of validity check

        rho(p.I) = rhoL;
        velx(p.I) = v_up(0);
        vely(p.I) = v_up(1);
        velz(p.I) = v_up(2);
        eps(p.I) = epsL;
        press(p.I) = pressL;
        entropy(p.I) = entropyL;
        Ye(p.I) = YeL;
        temperature(p.I) = tempL;

        saved_rho(p.I) = rhoL;
        saved_velx(p.I) = v_up(0);
        saved_vely(p.I) = v_up(1);
        saved_velz(p.I) = v_up(2);
        saved_eps(p.I) = epsL;
        saved_Ye(p.I) = YeL;

        CCTK_REAL wlor = calc_wlorentz(v_low, v_up);

        zvec_x(p.I) = wlor * v_up(0);
        zvec_y(p.I) = wlor * v_up(1);
        zvec_z(p.I) = wlor * v_up(2);

        svec_x(p.I) = (rhoL + rhoL * epsL + pressL) * wlor * wlor * v_up(0);
        svec_y(p.I) = (rhoL + rhoL * epsL + pressL) * wlor * wlor * v_up(1);
        svec_z(p.I) = (rhoL + rhoL * epsL + pressL) * wlor * wlor * v_up(2);
      });
}

extern "C" void AsterX_CheckPrims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CheckPrims;
  DECLARE_CCTK_PARAMETERS;

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
    // Get local eos object
    auto eos_1p_poly = global_eos_1p_poly;
    auto eos_3p_ig = global_eos_3p_ig;

    CheckPrims(cctkGH, eos_1p_poly, eos_3p_ig);
    break;
  }
  case eos_3param::Hybrid: {
    // Get local eos object
    auto eos_1p_poly = global_eos_1p_poly;
    auto eos_3p_hyb = global_eos_3p_hyb;

    CheckPrims(cctkGH, eos_1p_poly, eos_3p_hyb);
    break;
  }
  case eos_3param::Tabulated: {
    // Get local eos object
    auto eos_1p_poly = global_eos_1p_poly;
    auto eos_3p_tab3d = global_eos_3p_tab3d;

    CheckPrims(cctkGH, eos_1p_poly, eos_3p_tab3d);
    break;
  }
  default:
    assert(0);
  }
}

} // namespace AsterX
