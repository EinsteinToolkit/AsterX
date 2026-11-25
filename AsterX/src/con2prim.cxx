#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>

#include "c2p.hxx"
#include "c2p_1DEntropy.hxx"
#include "c2p_1DPalenzuela.hxx"
#include "c2p_1DRePrimAnd.hxx"
#include "c2p_2DNoble.hxx"

#include "aster_utils.hxx"
#include "atmo_global.hxx"
#include "setup_eos.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace EOSX;
using namespace Con2PrimFactory;
using namespace AsterUtils;

enum class eos_3param { IdealGas, Hybrid, Tabulated };
enum class c2p_first_t { None, Noble, RePrimAnd, Palenzuela, Entropy };
enum class c2p_second_t { None, Noble, RePrimAnd, Palenzuela, Entropy };

enum C2PFlag : CCTK_INT {
  C2P_INIT = 0,    // initial value
  C2P_PRIME = 1,   // first solver succeeded
  C2P_SECOND = 2,  // second solver succeeded
  C2P_ENTROPY = 3, // 1‑D Entropy (kappa) solver succeeded
  C2P_ATMO = 4,    // when (cv.dens <= sqrt_detg * rho_atmo_cut) is true
  C2P_AVG = 5,     // primitives obtained by neighbour‑averaging
  C2P_FAIL = 6     // when C2P fails
};

template <typename EOSIDType, typename EOSType>
void AsterX_Con2Prim_typeEoS(CCTK_ARGUMENTS, EOSIDType *eos_1p,
                             EOSType *eos_3p) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim;
  DECLARE_CCTK_PARAMETERS;

  const auto eos_1p_pwpoly = global_eos_1p_pwpoly;
  const void *eos_cold = eos_1p_pwpoly
                            ? static_cast<const void *>(eos_1p_pwpoly)
                            : eos_1p;
  atmosphere atmo_const{};
  const bool use_atmo_const =
      get_global_atmo(eos_cold, eos_3p, atmo_const);

  c2p_first_t c2p_fir;
  c2p_second_t c2p_sec;

  if (CCTK_EQUALS(c2p_prime, "Noble")) {
    c2p_fir = c2p_first_t::Noble;
  } else if (CCTK_EQUALS(c2p_prime, "RePrimAnd")) {
    c2p_fir = c2p_first_t::RePrimAnd;
  } else if (CCTK_EQUALS(c2p_prime, "Palenzuela")) {
    c2p_fir = c2p_first_t::Palenzuela;
  } else if (CCTK_EQUALS(c2p_prime, "Entropy")) {
    c2p_fir = c2p_first_t::Entropy;
  } else if (CCTK_EQUALS(c2p_prime, "None")) {
    c2p_fir = c2p_first_t::None;
  } else {
    CCTK_ERROR("Unknown value for parameter \"c2p_prime\"");
  }

  if (CCTK_EQUALS(c2p_second, "Noble")) {
    c2p_sec = c2p_second_t::Noble;
  } else if (CCTK_EQUALS(c2p_second, "RePrimAnd")) {
    c2p_sec = c2p_second_t::RePrimAnd;
  } else if (CCTK_EQUALS(c2p_second, "Palenzuela")) {
    c2p_sec = c2p_second_t::Palenzuela;
  } else if (CCTK_EQUALS(c2p_second, "Entropy")) {
    c2p_sec = c2p_second_t::Entropy;
  } else if (CCTK_EQUALS(c2p_second, "None")) {
    c2p_sec = c2p_second_t::None;
  } else {
    CCTK_ERROR("Unknown value for parameter \"c2p_second\"");
  }

  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  const auto c2p_impl = [=] CCTK_DEVICE(
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    // Note that HydroBaseX gfs are NaN when entering this loop due
    // explicit dependence on conservatives from
    // AsterX -> dependents tag

    // Setting up atmosphere
    atmosphere atmo{};
    if (use_atmo_const) {
      atmo = atmo_const;
    } else {
      const CCTK_REAL radial_distance = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
      if (eos_1p_pwpoly) {
        atmo = make_atmo(
            eos_1p_pwpoly, eos_3p, radial_distance, rho_abs_min, p_atmo,
            t_atmo, Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo,
            n_temp_atmo, atmo_tol, thermal_eos_atmo, use_press_atmo,
            Ye_atmo_beq);
      } else {
        atmo = make_atmo(
            eos_1p, eos_3p, radial_distance, rho_abs_min, p_atmo,
            t_atmo, Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo,
            n_temp_atmo, atmo_tol, thermal_eos_atmo, use_press_atmo,
            Ye_atmo_beq);
      }
    }
    const CCTK_REAL rho_atm = atmo.rho_atmo;
    const CCTK_REAL press_atm = atmo.press_atmo;
    const CCTK_REAL rho_atmo_cut = atmo.rho_cut;

    // ----- Construct C2P objects -----

    // Construct Noble c2p object:
    c2p_2DNoble c2p_Noble(eos_3p, atmo, max_iter, c2p_tol, alp_thresh, vw_lim,
                          B_lim, rho_BH, eps_BH, vwlim_BH, sigma_max,
                          inv_beta_max, Ye_lenient, use_z, use_temperature,
                          use_press_atmo, soft_root_convergence,
                          soft_root_width_factor);

    // Construct RePrimAnd c2p object:
    c2p_1DRePrimAnd c2p_RPA(eos_3p, atmo, max_iter, c2p_tol, alp_thresh, vw_lim,
                            B_lim, rho_BH, eps_BH, vwlim_BH, sigma_max,
                            inv_beta_max, Ye_lenient, use_z, use_temperature,
                            use_press_atmo, soft_root_convergence,
                            soft_root_width_factor);

    // Construct Palenzuela c2p object:
    c2p_1DPalenzuela c2p_Pal(eos_3p, atmo, max_iter, c2p_tol, alp_thresh,
                             vw_lim, B_lim, rho_BH, eps_BH, vwlim_BH, sigma_max,
                             inv_beta_max, Ye_lenient, use_z, use_temperature,
                             use_press_atmo, soft_root_convergence,
                             soft_root_width_factor);

    // Construct Entropy c2p object:
    c2p_1DEntropy c2p_Ent(eos_3p, atmo, max_iter, c2p_tol, alp_thresh, vw_lim,
                          B_lim, rho_BH, eps_BH, vwlim_BH, sigma_max,
                          inv_beta_max, Ye_lenient, use_z, use_temperature,
                          use_press_atmo, soft_root_convergence,
                          soft_root_width_factor);

    // ----------

    /* Get lapse */
    const CCTK_REAL alp_avg = calc_avg_v2c(alp, p);

    /* Get shift */
    const vec<CCTK_REAL, 3> beta_avg(
        [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p); });

    /* Get covariant metric */
    const smat<CCTK_REAL, 3> glo(
        [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });

    /* Get mask */
    CCTK_REAL mask_local = use_mask ? calc_avg_v2c(aster_mask_vc, p) : 1.0;

    /* Calculate inverse of 3-metric */
    const CCTK_REAL spatial_detg = calc_det(glo);
    const CCTK_REAL sqrt_detg = sqrt(spatial_detg);

    vec<CCTK_REAL, 3> v_up{saved_velx(p.I), saved_vely(p.I), saved_velz(p.I)};
    vec<CCTK_REAL, 3> v_low = calc_contraction(glo, v_up);
    CCTK_REAL zsq{0.0};

    // TODO: Debug code to capture v>1 early,
    // remove soon
    const CCTK_REAL vsq = calc_contraction(v_low, v_up);
    if (vsq >= 1.0) {
      CCTK_REAL wlim = sqrt(1.0 + vw_lim * vw_lim);
      CCTK_REAL vlim = vw_lim / wlim;
      v_up *= vlim / sqrt(vsq);
      v_low *= vlim / sqrt(vsq);
      zsq = vw_lim * vw_lim;
    } else {
      zsq = vsq / (1.0 - vsq);
    }

    // CCTK_REAL wlor = calc_wlorentz(v_low, v_up);
    CCTK_REAL wlor = sqrt(1.0 + zsq);

    // Note that cv are densitized, i.e. they all include sqrt_detg
    cons_vars cv{dens(p.I), {momx(p.I), momy(p.I), momz(p.I)},
                 tau(p.I),  DYe(p.I),
                 DEnt(p.I), {dBx(p.I), dBy(p.I), dBz(p.I)}};

    // Undensitized magnetic fields
    const vec<CCTK_REAL, 3> Bup{cv.dBvec(0) / sqrt_detg,
                                cv.dBvec(1) / sqrt_detg,
                                cv.dBvec(2) / sqrt_detg};

    prim_vars pv;
    prim_vars pv_seeds{saved_rho(p.I),
                       saved_eps(p.I),
                       saved_Ye(p.I),
                       eos_3p->press_from_rho_eps_ye(
                           saved_rho(p.I), saved_eps(p.I), saved_Ye(p.I)),
                       temperature(p.I),
                       eos_3p->kappa_from_rho_eps_ye(
                           saved_rho(p.I), saved_eps(p.I), saved_Ye(p.I)),
                       v_up,
                       wlor,
                       Bup};

    /* set flag to success */
    bool c2p_flag_local = true;
    CCTK_INT c2p_flag_code = C2P_INIT;
    bool call_c2p = true;
    bool bh_failed = false;

    // Check if point is below atmosphere, and if atmosphere obeys magnetization
    // limits (RPA only). Magnetization limits are currently only applied for RPA C2P, 
    // while they are not obeyed in the other cases in the atmopshere -> TODO
    const CCTK_REAL b2_atm = calc_norm(Bup, glo);
    const bool set_atmo = (cv.dens <= sqrt_detg * rho_atmo_cut) &&
                          (c2p_off_floor_strict ||
                           ((b2_atm / rho_atm <= sigma_max) &&
                            (b2_atm / (2 * press_atm) <= inv_beta_max)));
    if (set_atmo) {
      pv.Bvec = Bup;
      atmo.set(pv, cv, glo);
      atmo.set(pv_seeds);
      c2p_flag_code = C2P_ATMO;
      call_c2p = false;
    }

    // Modifying primitive seeds within BH interiors before C2Ps are called
    // NOTE: By default, alp_thresh=0 so the if condition below is never
    // triggered. One must be very careful when using this functionality and
    // must correctly set alp_thresh, rho_BH, eps_BH and vwlim_BH in the parfile

    if (alp_avg < alp_thresh) {
      mask_local = 0.0;
    }
    aster_mask_cc(p.I) = mask_local;

    if (excise) {

      if (mask_local != 1.0) {
        bh_failed =
            !c2p_Noble.bh_interior<EOSType, false>(eos_3p, pv_seeds, cv, glo);
        pv = pv_seeds;
        call_c2p = false;
      }
    }

    // Construct error report object:
    c2p_report rep_first;
    c2p_report rep_second;
    c2p_report rep_ent;

    // ----- ----- C2P ----- -----

    if (call_c2p) {
      // Limit conservatives before calling C2P
      // Do not modify a completed atmosphere or excision reset.
      c2p_Noble.cons_floors_and_ceilings(eos_3p, cv, glo, tauFluid_atmo);

      // Calling the first C2P
      c2p_flag_code = C2P_PRIME;
      switch (c2p_fir) {
      case c2p_first_t::Noble: {
        c2p_Noble.solve(eos_3p, pv, pv_seeds, cv, alp_avg, beta_avg, glo,
                        rep_first, use_entropy_fix);
        break;
      }
      case c2p_first_t::RePrimAnd: {
        c2p_RPA.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_first,
                      use_entropy_fix);
        break;
      }
      case c2p_first_t::Palenzuela: {
        c2p_Pal.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_first,
                      use_entropy_fix);
        break;
      }
      case c2p_first_t::Entropy: {
        c2p_Ent.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_first);
        break;
      }
      case c2p_first_t::None: {
        // solve not called, pv remains unwritten
        break;
      }
      default:
        assert(0);
      }

      if (rep_first.failed()) {
        c2p_flag_code = C2P_SECOND;

        if (debug_mode) {
          printf("First C2P failed :( \n");
          rep_first.debug_message();
          printf("Calling the back up C2P.. \n");
        }

        // Calling the second C2P
        switch (c2p_sec) {
        case c2p_second_t::Noble: {
          c2p_Noble.solve(eos_3p, pv, pv_seeds, cv, alp_avg, beta_avg, glo,
                          rep_second, use_entropy_fix);
          break;
        }
        case c2p_second_t::RePrimAnd: {
          c2p_RPA.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_second,
                        use_entropy_fix);
          break;
        }
        case c2p_second_t::Palenzuela: {
          c2p_Pal.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_second,
                        use_entropy_fix);
          break;
        }
        case c2p_second_t::Entropy: {
          c2p_Ent.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_second);
          break;
        }
        case c2p_second_t::None: {
          // solve not called, pv remains unwritten
          break;
        }
        default:
          assert(0);
        }
      }

      if (rep_first.failed() && rep_second.failed()) {

        if (use_entropy_fix) {

          c2p_flag_code = C2P_ENTROPY;
          c2p_Ent.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_ent);

          if (rep_ent.failed()) {

            c2p_flag_local = false;
            c2p_flag_code = C2P_FAIL;

            if (debug_mode) {
              printf("Entropy C2P failed. Setting point to atmosphere.\n");
              rep_ent.debug_message();
              printf("WARNING: \n"
                     "C2Ps failed. Printing cons and saved prims before set to "
                     "atmo: \n"
                     "cctk_iteration = %i \n "
                     "x, y, z = %26.16e, %26.16e, %26.16e \n "
                     "dens = %26.16e \n tau = %26.16e \n momx = %26.16e \n "
                     "momy = %26.16e \n momz = %26.16e \n dBx = %26.16e \n "
                     "dBy = %26.16e \n dBz = %26.16e \n "
                     "saved_rho = %26.16e \n saved_eps = %26.16e \n press= "
                     "%26.16e "
                     "\n "
                     "saved_velx = %26.16e \n saved_vely = %26.16e \n "
                     "saved_velz = "
                     "%26.16e \n "
                     "Bvecx = %26.16e \n Bvecy = %26.16e \n "
                     "Bvecz = %26.16e \n "
                     "Avec_x = %26.16e \n Avec_y = %26.16e \n Avec_z = %26.16e "
                     "\n ",
                     cctk_iteration, p.x, p.y, p.z, dens(p.I), tau(p.I),
                     momx(p.I), momy(p.I), momz(p.I), dBx(p.I), dBy(p.I),
                     dBz(p.I), pv.rho, pv.eps, pv.press, pv.vel(0), pv.vel(1),
                     pv.vel(2), pv.Bvec(0), pv.Bvec(1), pv.Bvec(2),
                     // rho(p.I), eps(p.I), press(p.I), velx(p.I), vely(p.I),
                     // velz(p.I), Bvecx(p.I), Bvecy(p.I), Bvecz(p.I),
                     Avec_x(p.I), Avec_y(p.I), Avec_z(p.I));
            }

            if (mask_local != 1.0) {
              // Failure inside mask
              bh_failed =
                  !c2p_Noble.bh_interior<EOSType, false>(eos_3p, pv_seeds, cv, glo);
              pv = pv_seeds;
            } else {
              // Failure outside, set to atmo
              cv.dBvec(0) = sqrt_detg * Bup(0);
              cv.dBvec(1) = sqrt_detg * Bup(1);
              cv.dBvec(2) = sqrt_detg * Bup(2);
              pv.Bvec = Bup;
              atmo.set(pv, cv, glo);
            }
          }

        } else {

          c2p_flag_local = false;
          c2p_flag_code = C2P_FAIL;

          if (debug_mode) {
            printf("Second C2P failed too :( :( \n");
            rep_second.debug_message();
            printf(
                "WARNING: \n"
                "C2Ps failed. Printing cons and saved prims before set to "
                "atmo: \n"
                "cctk_iteration = %i \n "
                "x, y, z = %26.16e, %26.16e, %26.16e \n "
                "dens = %26.16e \n tau = %26.16e \n momx = %26.16e \n "
                "momy = %26.16e \n momz = %26.16e \n dBx = %26.16e \n "
                "dBy = %26.16e \n dBz = %26.16e \n "
                "saved_rho = %26.16e \n saved_eps = %26.16e \n press= %26.16e "
                "\n "
                "saved_velx = %26.16e \n saved_vely = %26.16e \n saved_velz = "
                "%26.16e \n "
                "Bvecx = %26.16e \n Bvecy = %26.16e \n "
                "Bvecz = %26.16e \n "
                "Avec_x = %26.16e \n Avec_y = %26.16e \n Avec_z = %26.16e \n ",
                cctk_iteration, p.x, p.y, p.z, dens(p.I), tau(p.I), momx(p.I),
                momy(p.I), momz(p.I), dBx(p.I), dBy(p.I), dBz(p.I), pv.rho,
                pv.eps, pv.press, pv.vel(0), pv.vel(1), pv.vel(2), pv.Bvec(0),
                pv.Bvec(1), pv.Bvec(2),
                // rho(p.I), eps(p.I), press(p.I), velx(p.I), vely(p.I),
                // velz(p.I), Bvecx(p.I), Bvecy(p.I), Bvecz(p.I),
                Avec_x(p.I), Avec_y(p.I), Avec_z(p.I));
          }

          if (mask_local != 1.0) {
            // Failure inside mask
            bh_failed =
                !c2p_Noble.bh_interior<EOSType, false>(eos_3p, pv_seeds, cv, glo);
            pv = pv_seeds;
          } else {
            // Failure outside, set to atmo
            cv.dBvec(0) = sqrt_detg * Bup(0);
            cv.dBvec(1) = sqrt_detg * Bup(1);
            cv.dBvec(2) = sqrt_detg * Bup(2);
            pv.Bvec = Bup;
            atmo.set(pv, cv, glo);
          }
        }
      }

      // Inside mask, C2P success
      if ((mask_local != 1.0) && c2p_flag_local) {
        bh_failed = !c2p_Noble.bh_interior<EOSType, true>(eos_3p, pv, cv, glo);
      }
    }

    if (bh_failed) {
      pv.Bvec = Bup;
      cv.dBvec = sqrt_detg * Bup;
      atmo.set(pv, cv, glo);
      c2p_flag_code = C2P_FAIL;
    }

    con2prim_flag(p.I) = c2p_flag_code;

    // ----- ----- C2P ----- -----

    // ----- Write to gfs -----

    // dummy vars
    CCTK_REAL Ex, Ey, Ez;

    // Write back pv
    pv.scatter(rho(p.I), eps(p.I), Ye(p.I), press(p.I), temperature(p.I),
               entropy(p.I), velx(p.I), vely(p.I), velz(p.I), wlor, Bvecx(p.I),
               Bvecy(p.I), Bvecz(p.I), Ex, Ey, Ez);

    zvec_x(p.I) = wlor * pv.vel(0);
    zvec_y(p.I) = wlor * pv.vel(1);
    zvec_z(p.I) = wlor * pv.vel(2);

    svec_x(p.I) =
        (pv.rho + pv.rho * pv.eps + pv.press) * wlor * wlor * pv.vel(0);
    svec_y(p.I) =
        (pv.rho + pv.rho * pv.eps + pv.press) * wlor * wlor * pv.vel(1);
    svec_z(p.I) =
        (pv.rho + pv.rho * pv.eps + pv.press) * wlor * wlor * pv.vel(2);

    // Write back cv
    cv.scatter(dens(p.I), momx(p.I), momy(p.I), momz(p.I), tau(p.I), DYe(p.I),
               DEnt(p.I), dBx(p.I), dBy(p.I), dBz(p.I));

    // Update saved prims
    saved_rho(p.I) = rho(p.I);
    saved_velx(p.I) = velx(p.I);
    saved_vely(p.I) = vely(p.I);
    saved_velz(p.I) = velz(p.I);
    saved_eps(p.I) = eps(p.I);
    saved_Ye(p.I) = Ye(p.I);

    // Compute diagnostics
    v_low = calc_contraction(glo, pv.vel);
    const vec<CCTK_REAL, 3> B_low = calc_contraction(glo, pv.Bvec);

    const CCTK_REAL alp_b0 = wlor * calc_contraction(pv.Bvec, v_low);

    const CCTK_REAL B2 = calc_contraction(pv.Bvec, B_low);
    const CCTK_REAL bsq = (B2 + alp_b0 * alp_b0) / (wlor * wlor);

    w_lorentz(p.I) = wlor;
    B_norm(p.I) = sqrt(B2);
    b2small(p.I) = bsq;
    volform(p.I) = sqrt_detg;
  };

  cctk_grid.loop_all_device<1, 1, 1>(grid.nghostzones, c2p_impl);
}

extern "C" void AsterX_Con2Prim(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Con2Prim;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(evolution_eos, "Hybrid") && thermal_eos_atmo) {
    CCTK_ERROR("Hybrid EOS does not implement *_from_rho_temp_ye; set "
               "Con2PrimFactory::thermal_eos_atmo = no.");
  }
  if (CCTK_EQUALS(evolution_eos, "Tabulated3d") && !use_temperature) {
    CCTK_ERROR("Tabulated3d requires Con2PrimFactory::use_temperature = yes.");
  }
  if (CCTK_EQUALS(evolution_eos, "Tabulated3d") && use_press_atmo) {
    CCTK_ERROR("Tabulated3d does not support eps_from_rho_press_ye; set "
               "Con2PrimFactory::use_press_atmo = no.");
  }

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
    // Get local eos objects
    auto eos_1p_poly = global_eos_1p_poly;
    auto eos_3p_ig = global_eos_3p_ig;

    AsterX_Con2Prim_typeEoS(CCTK_PASS_CTOC, eos_1p_poly, eos_3p_ig);
    break;
  }
  case eos_3param::Hybrid: {
    if (global_eos_3p_hyb_pwpoly) {
      // pwpoly cold + pwpoly-hybrid
      auto eos_cold = global_eos_1p_pwpoly;
      auto eos_3p_hyb = global_eos_3p_hyb_pwpoly;

      if (!eos_cold) {
        CCTK_ERROR("Hybrid(PWPolytropic) selected but no pwpoly cold EOS was "
                   "initialized");
      }
      AsterX_Con2Prim_typeEoS(CCTK_PASS_CTOC, eos_cold, eos_3p_hyb);

    } else if (global_eos_3p_hyb_poly) {
      // poly cold + poly-hybrid
      auto eos_cold = global_eos_1p_poly;
      auto eos_3p_hyb = global_eos_3p_hyb_poly;

      if (!eos_cold) {
        CCTK_ERROR("Hybrid(Polytropic) selected but no polytropic cold EOS was "
                   "initialized");
      }
      AsterX_Con2Prim_typeEoS(CCTK_PASS_CTOC, eos_cold, eos_3p_hyb);

    } else {
      CCTK_ERROR(
          "Hybrid EOS selected but no hybrid EOS object was initialized");
    }

    break;
  }
  case eos_3param::Tabulated: {
    // Get local eos objects
    auto eos_1p_poly = global_eos_1p_poly;
    auto eos_3p_tab3d = global_eos_3p_tab3d;

    AsterX_Con2Prim_typeEoS(CCTK_PASS_CTOC, eos_1p_poly, eos_3p_tab3d);
    break;
  }
  default:
    assert(0);
  }
}

template <typename EOSIDType, typename EOSType>
void InterpolateFailed(CCTK_ARGUMENTS, const EOSIDType *eos_1p,
                       const EOSType *eos_3p) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim_Interpolate_Failed;
  DECLARE_CCTK_PARAMETERS;
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  grid.loop_int_device<1, 1, 1>(grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) {
        if (con2prim_flag(p.I) != CCTK_REAL(C2P_FAIL))
          return;
        // In-place updates can race between failed neighbours.
        const auto flag_nbs = get_neighbors(con2prim_flag, p);
        const auto rho_nbs = get_neighbors(rho, p);
        const auto eps_nbs = get_neighbors(eps, p);
        const auto Ye_nbs = get_neighbors(Ye, p);
        const auto velx_nbs = get_neighbors(velx, p);
        const auto vely_nbs = get_neighbors(vely, p);
        const auto velz_nbs = get_neighbors(velz, p);
        const auto good_nb = [&](int i) {
          return flag_nbs(i) != CCTK_REAL(C2P_FAIL) &&
                 flag_nbs(i) != CCTK_REAL(C2P_INIT) &&
                 std::isfinite(rho_nbs(i)) && std::isfinite(eps_nbs(i)) &&
                 std::isfinite(Ye_nbs(i)) && std::isfinite(velx_nbs(i)) &&
                 std::isfinite(vely_nbs(i)) && std::isfinite(velz_nbs(i));
        };
        CCTK_REAL sum_nbs = 0.0;
        for (int i = 0; i < 6; ++i)
          sum_nbs += good_nb(i);
        if (sum_nbs == 0.0)
          return;
        const auto average = [&](const auto &values) {
          CCTK_REAL value = 0.0;
          for (int i = 0; i < 6; ++i)
            if (good_nb(i))
              value += values(i);
          return value / sum_nbs;
        };
        const auto atmo = make_atmo(eos_1p, eos_3p,
            sqrt(p.x * p.x + p.y * p.y + p.z * p.z), rho_abs_min, p_atmo, t_atmo,
            Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo, n_temp_atmo, atmo_tol,
            thermal_eos_atmo, use_press_atmo, Ye_atmo_beq);
        const smat<CCTK_REAL, 3> g([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });
        const CCTK_REAL detg = calc_det(g);
        if (!std::isfinite(detg) || detg <= 0.0 || g(0, 0) <= 0.0 ||
            g(0, 0) * g(1, 1) <= g(0, 1) * g(0, 1))
          return;
        prim_vars pv;
        pv.Bvec = {Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};
        if (!std::isfinite(pv.Bvec(0)) || !std::isfinite(pv.Bvec(1)) ||
            !std::isfinite(pv.Bvec(2)))
          return;
        const CCTK_REAL rhoL = average(rho_nbs);
        if (rhoL <= atmo.rho_cut) {
          atmo.set(pv);
        } else {
          pv.rho = rhoL;
          pv.eps = average(eps_nbs);
          pv.Ye = average(Ye_nbs);
          if (!complete_prims(eos_3p, pv, false))
            return;
          if (use_press_atmo && pv.press < atmo.press_atmo) {
            pv.eps = eos_3p->eps_from_rho_press_ye(pv.rho, atmo.press_atmo, pv.Ye);
            if (!complete_prims(eos_3p, pv, false))
              return;
          } else if (!use_press_atmo && pv.temperature < atmo.temp_atmo) {
            pv.temperature = atmo.temp_atmo;
            if (!complete_prims(eos_3p, pv, true))
              return;
          }
          pv.vel = {average(velx_nbs), average(vely_nbs), average(velz_nbs)};
          const auto v_low = calc_contraction(g, pv.vel);
          const CCTK_REAL vsq = calc_contraction(v_low, pv.vel);
          if (!std::isfinite(vsq) || vsq < 0.0)
            return;
          const CCTK_REAL vlim = vw_lim / sqrt(1.0 + vw_lim * vw_lim);
          if (vsq > vlim * vlim)
            pv.vel *= vlim / sqrt(vsq);
          pv.w_lor = calc_wlorentz(calc_contraction(g, pv.vel), pv.vel);
          pv.E = calc_contraction(calc_inv(g, detg),
                                   calc_cross_product(pv.Bvec, pv.vel));
        }
        cons_vars cv;
        cv.from_prim(pv, g);
        if (!std::isfinite(cv.dens) || cv.dens <= 0.0 ||
            !std::isfinite(cv.tau) || !std::isfinite(cv.DYe) ||
            !std::isfinite(cv.DEnt))
          return;
        for (int d = 0; d < 3; ++d)
          if (!std::isfinite(cv.mom(d)) || !std::isfinite(cv.dBvec(d)))
            return;
        CCTK_REAL Ex, Ey, Ez, Bx, By, Bz;
        pv.scatter(rho(p.I), eps(p.I), Ye(p.I), press(p.I), temperature(p.I),
            entropy(p.I), velx(p.I), vely(p.I), velz(p.I), w_lorentz(p.I),
            Bx, By, Bz, Ex, Ey, Ez);
        // The magnetic field is unchanged, so do not write the staggered B.
        dens(p.I) = cv.dens;
        momx(p.I) = cv.mom(0);
        momy(p.I) = cv.mom(1);
        momz(p.I) = cv.mom(2);
        tau(p.I) = cv.tau;
        DYe(p.I) = cv.DYe;
        DEnt(p.I) = cv.DEnt;
        saved_rho(p.I) = pv.rho;
        saved_eps(p.I) = pv.eps;
        saved_Ye(p.I) = pv.Ye;
        saved_velx(p.I) = pv.vel(0);
        saved_vely(p.I) = pv.vel(1);
        saved_velz(p.I) = pv.vel(2);
        const CCTK_REAL rhoh = pv.rho * (1.0 + pv.eps) + pv.press;
        zvec_x(p.I) = pv.w_lor * pv.vel(0);
        zvec_y(p.I) = pv.w_lor * pv.vel(1);
        zvec_z(p.I) = pv.w_lor * pv.vel(2);
        svec_x(p.I) = rhoh * pv.w_lor * pv.w_lor * pv.vel(0);
        svec_y(p.I) = rhoh * pv.w_lor * pv.w_lor * pv.vel(1);
        svec_z(p.I) = rhoh * pv.w_lor * pv.w_lor * pv.vel(2);
        const CCTK_REAL B2 = calc_contraction(pv.Bvec, calc_contraction(g, pv.Bvec));
        const CCTK_REAL Bv = calc_contraction(pv.Bvec, calc_contraction(g, pv.vel));
        B_norm(p.I) = sqrt(B2);
        b2small(p.I) = B2 / (pv.w_lor * pv.w_lor) + Bv * Bv;
        volform(p.I) = sqrt(detg);
        con2prim_flag(p.I) = C2P_AVG;
      });
}

extern "C" void AsterX_Con2Prim_Interpolate_Failed(CCTK_ARGUMENTS) {
  const auto run = [&](const auto *eos) {
    if (global_eos_1p_pwpoly)
      InterpolateFailed(cctkGH, global_eos_1p_pwpoly, eos);
    else
      InterpolateFailed(cctkGH, global_eos_1p_poly, eos);
  };
  if (global_eos_3p_ig)
    run(global_eos_3p_ig);
  else if (global_eos_3p_tab3d)
    run(global_eos_3p_tab3d);
  else
    CCTK_ERROR("Neighbour C2P repair supports IdealGas and Tabulated3d");
}

} // namespace AsterX
