#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

#include "aster_utils.hxx"

#include "reprimand/eos_thermal.h" // The EOS framework
#include "reprimand/eos_idealgas.h"
#include "reprimand/eos_barotropic.h"
#include "reprimand/eos_barotr_poly.h"
#include "reprimand/con2prim_imhd.h" // The con2prim framework

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace EOS_Toolkit;
using namespace AsterUtils;

extern "C" void AsterX_Con2Prim(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim;
  DECLARE_CCTK_PARAMETERS;

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  // Loop over the interior of the grid
  cctk_grid.loop_int_device<
      1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                     const PointDesc
                                         &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    // Setting up initial data EOS
    const CCTK_REAL n = 1 / (poly_gamma - 1); // Polytropic index
    const CCTK_REAL adiab_ind_id = 1.0 / (poly_gamma - 1);
    const CCTK_REAL rmd_p = pow(poly_k, -n); // Polytropic density scale
    const auto eos_id = make_eos_barotr_poly(adiab_ind_id, rmd_p, rho_max);

    // Setting up evolution EOS
    const CCTK_REAL adiab_ind_evol = 1.0 / (gl_gamma - 1);
    const auto eos = make_eos_idealgas(adiab_ind_evol, eps_max, rho_max);

    // Setting up atmosphere

    CCTK_REAL rho_atm = 0.0;   // dummy initialization
    CCTK_REAL press_atm = 0.0; // dummy initialization
    CCTK_REAL eps_atm = 0.0;   // dummy initialization
    CCTK_REAL radial_distance = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);

    // Grading rho
    rho_atm = (radial_distance > r_atmo)
                  ? (rho_abs_min * pow((r_atmo / radial_distance), n_rho_atmo))
                  : rho_abs_min;
    const CCTK_REAL rho_atmo_cut = rho_atm * (1 + atmo_tol);

    // Grading pressure and eps
    if (thermal_eos_atmo) {
      press_atm = (radial_distance > r_atmo)
                      ? (p_atmo * pow(r_atmo / radial_distance, n_press_atmo))
                      : p_atmo;
      // TODO: eos.at_rho_press_ye(rho_atm, press_atm, Ye_atmo).eps() does not
      // exist in RePrimAnd Currently computing eps from ideal gas EOS eps_atm =
      // eos.at_rho_press_ye(rho_atm, press_atm, Ye_atmo).eps();
      eps_atm = press_atm / (rho_atm * (gl_gamma - 1.));
    } else {
      eps_atm = eos_id.at_rho(rho_atm).eps();
      eps_atm = eos.range_eps(rho_atm, Ye_atmo).limit_to(eps_atm);
      press_atm = eos.at_rho_eps_ye(rho_atm, eps_atm, Ye_atmo).press();
    }
    const atmosphere atmo(rho_atm, eps_atm, Ye_atmo, press_atm, rho_atmo_cut);

    CCTK_REAL dummy_Ye = 0.5;
    CCTK_REAL dummy_dYe = 0.5;

    // Get a recovery function
    con2prim_mhd cv2pv(eos, rho_strict, Ye_lenient, vw_lim, B_lim, atmo,
                       c2p_tol, max_iter);

    /* Get covariant metric */
    const smat<CCTK_REAL, 3> glo(
        [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });

    sm_metric3 g(sm_symt3l(glo(0, 0), glo(0, 1), glo(1, 1), glo(0, 2),
                           glo(1, 2), glo(2, 2)));

    /* Get covariant metric used for BH interior fixes */
    const smat<CCTK_REAL, 3> glow(
        [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });

    /* Calculate inverse of 3-metric */
    const CCTK_REAL spatial_detg = calc_det(glo);
    const CCTK_REAL sqrt_detg = sqrt(spatial_detg);

    prim_vars_mhd pv;
    prim_vars_mhd pv_seeds{saved_rho(p.I),
                           saved_eps(p.I),
                           dummy_Ye,
                           press(p.I),
                           {saved_velx(p.I), saved_vely(p.I), saved_velz(p.I)},
                           dummy_Ye,
                           {dBx(p.I), dBy(p.I), dBz(p.I)},
                           {dBx(p.I), dBy(p.I), dBz(p.I)}};
    con2prim_mhd::report rep;
    // Note that cv are densitized, i.e. they all include sqrt_detg
    cons_vars_mhd cv{dens(p.I),
                     tau(p.I),
                     dummy_dYe,
                     {momx(p.I), momy(p.I), momz(p.I)},
                     {dBx(p.I), dBy(p.I), dBz(p.I)}};

    // Modifying primitive seeds within BH interiors before C2Ps are called
    // NOTE: By default, alp_thresh=0 so the if condition below is never
    // triggered. One must be very careful when using this functionality and
    // must correctly set alp_thresh, rho_BH, eps_BH and vwlim_BH in the parfile

    if (alp(p.I) < alp_thresh) {
      if ((pv_seeds.rho > rho_BH) || (pv_seeds.eps > eps_BH)) {
        pv_seeds.rho = rho_BH; // typically set to 0.01% to 1% of rho_max of
                               // initial NS or disk
        pv_seeds.eps = eps_BH;
        pv_seeds.ye = Ye_atmo;
        pv_seeds.press = eos.at_rho_eps_ye(rho_BH, eps_BH, Ye_atmo).press();
        // check on velocities
        CCTK_REAL wlim_BH = sqrt(1.0 + vwlim_BH * vwlim_BH);
        CCTK_REAL vlim_BH = vwlim_BH / wlim_BH;
        CCTK_REAL sol_v =
            sqrt((pv_seeds.w_lor * pv_seeds.w_lor - 1.0)) / pv_seeds.w_lor;
        if (sol_v > vlim_BH) {
          pv_seeds.vel *= vlim_BH / sol_v;
          pv_seeds.w_lor = wlim_BH;
        }

        // cv.from_prim(pv_seeds, g);
        // We do not save electric field, required by cv.from_prim(pv_seeds, g)
        // in RPA. Thus, we recompute the CVs explicitly below
        const vec<CCTK_REAL, 3> &v_up = pv_seeds.vel;
        const vec<CCTK_REAL, 3> v_low = calc_contraction(glow, v_up);
        /* Computing B_j */
        const vec<CCTK_REAL, 3> &B_up = pv_seeds.B;
        const vec<CCTK_REAL, 3> B_low = calc_contraction(glow, B_up);
        /* Computing b^t : this is b^0 * alp */
        const CCTK_REAL bst = pv_seeds.w_lor * calc_contraction(B_up, v_low);
        /* Computing b_j */
        const vec<CCTK_REAL, 3> b_low = B_low / pv_seeds.w_lor + bst * v_low;
        /* Computing b^mu b_mu */
        const CCTK_REAL bs2 = (calc_contraction(B_up, B_low) + bst * bst) /
                              (pv_seeds.w_lor * pv_seeds.w_lor);
        // computing conserved from primitives
        cv.dens = sqrt_detg * pv_seeds.rho * pv_seeds.w_lor;
        cv.scon(0) =
            sqrt_detg *
            (pv_seeds.w_lor * pv_seeds.w_lor *
                 (pv_seeds.rho * (1.0 + pv_seeds.eps) + pv_seeds.press + bs2) *
                 v_low(0) -
             bst * b_low(0));
        cv.scon(1) =
            sqrt_detg *
            (pv_seeds.w_lor * pv_seeds.w_lor *
                 (pv_seeds.rho * (1.0 + pv_seeds.eps) + pv_seeds.press + bs2) *
                 v_low(1) -
             bst * b_low(1));
        cv.scon(2) =
            sqrt_detg *
            (pv_seeds.w_lor * pv_seeds.w_lor *
                 (pv_seeds.rho * (1.0 + pv_seeds.eps) + pv_seeds.press + bs2) *
                 v_low(2) -
             bst * b_low(2));
        cv.tau = sqrt_detg * (pv_seeds.w_lor * pv_seeds.w_lor *
                                  (pv_seeds.rho * (1.0 + pv_seeds.eps) +
                                   pv_seeds.press + bs2) -
                              (pv_seeds.press + 0.5 * bs2) - bst * bst) -
                 cv.dens;
        cv.bcons = sqrt_detg * pv_seeds.B;
        cv.tracer_ye = cv.dens * pv_seeds.ye;
      }
    }

    // Calling C2P
    cv2pv(pv, cv, g, rep);

    // Handle incorrectable errors
    if (rep.failed()) {
      CCTK_WARN(1, rep.debug_message().c_str());

      if (debug_mode) {
        // need to fix pv to computed values like pv.rho instead of rho(p.I)
        printf(
            "WARNING: "
            "C2P failed. Printing cons and saved prims before set to atmo or "
            "BH interior fix: \n"
            "cctk_iteration = %i \n "
            "x, y, z = %26.16e, %26.16e, %26.16e \n "
            "gxx, gxy, gxz, gyy, gyz, gzz = %f, %f, %f, %f, %f, %f \n "
            "dens = %26.16e \n tau = %26.16e \n momx = %26.16e \n "
            "momy = %26.16e \n momz = %26.16e \n dBx = %26.16e \n "
            "dBy = %26.16e \n dBz = %26.16e \n "
            "saved_rho = %26.16e \n saved_eps = %26.16e \n press= %26.16e \n "
            "saved_velx = %26.16e \n saved_vely = %26.16e \n saved_velz = "
            "%26.16e \n "
            "Bvecx = %26.16e \n Bvecy = %26.16e \n "
            "Bvecz = %26.16e \n "
            "Avec_x = %26.16e \n Avec_y = %26.16e \n Avec_z = %26.16e \n ",
            cctk_iteration, p.x, p.y, p.z, glo(0, 0), glo(0, 1), glo(0, 2),
            glo(1, 1), glo(1, 2), glo(2, 2), dens(p.I), tau(p.I), momx(p.I),
            momy(p.I), momz(p.I), dBx(p.I), dBy(p.I), dBz(p.I), pv.rho, pv.eps,
            pv.press, pv.vel(0), pv.vel(1), pv.vel(2), pv.B(0), pv.B(1),
            pv.B(2),
            // rho(p.I), eps(p.I), press(p.I), velx(p.I), vely(p.I),
            // velz(p.I), Bvecx(p.I), Bvecy(p.I), Bvecz(p.I),
            Avec_x(p.I), Avec_y(p.I), Avec_z(p.I));
      }

      if (alp(p.I) < alp_thresh) {
        if ((pv.rho > rho_BH) || (pv.eps > eps_BH)) {
          pv.rho = rho_BH; // typically set to 0.01% to 1% of rho_max of
                           // initial NS or disk
          pv.eps = eps_BH;
          pv.ye = Ye_atmo;
          pv.press = eos.at_rho_eps_ye(rho_BH, eps_BH, Ye_atmo).press();
          // check on velocities
          CCTK_REAL wlim_BH = sqrt(1.0 + vwlim_BH * vwlim_BH);
          CCTK_REAL vlim_BH = vwlim_BH / wlim_BH;
          CCTK_REAL sol_v = sqrt((pv.w_lor * pv.w_lor - 1.0)) / pv.w_lor;
          if (sol_v > vlim_BH) {
            pv.vel *= vlim_BH / sol_v;
            pv.w_lor = wlim_BH;
          }
          // cv.from_prim(pv, g);
          // We do not save electric field, required by cv.from_prim(pv, g) in
          // RPA. Thus, we recompute the CVs explicitly below

          const vec<CCTK_REAL, 3> &v_up = pv.vel;
          const vec<CCTK_REAL, 3> v_low = calc_contraction(glow, v_up);
          /* Computing B_j */
          const vec<CCTK_REAL, 3> &B_up = pv.B;
          const vec<CCTK_REAL, 3> B_low = calc_contraction(glow, B_up);
          /* Computing b^t : this is b^0 * alp */
          const CCTK_REAL bst = pv.w_lor * calc_contraction(B_up, v_low);
          /* Computing b_j */
          const vec<CCTK_REAL, 3> b_low = B_low / pv.w_lor + bst * v_low;
          /* Computing b^mu b_mu */
          const CCTK_REAL bs2 = (calc_contraction(B_up, B_low) + bst * bst) /
                                (pv.w_lor * pv.w_lor);
          // computing conserved from primitives
          cv.dens = sqrt_detg * pv.rho * pv.w_lor;
          cv.scon(0) =
              sqrt_detg *
              (pv.w_lor * pv.w_lor *
                   (pv.rho * (1.0 + pv.eps) + pv.press + bs2) * v_low(0) -
               bst * b_low(0));
          cv.scon(1) =
              sqrt_detg *
              (pv.w_lor * pv.w_lor *
                   (pv.rho * (1.0 + pv.eps) + pv.press + bs2) * v_low(1) -
               bst * b_low(1));
          cv.scon(2) =
              sqrt_detg *
              (pv.w_lor * pv.w_lor *
                   (pv.rho * (1.0 + pv.eps) + pv.press + bs2) * v_low(2) -
               bst * b_low(2));
          cv.tau = sqrt_detg * (pv.w_lor * pv.w_lor *
                                    (pv.rho * (1.0 + pv.eps) + pv.press + bs2) -
                                (pv.press + 0.5 * bs2) - bst * bst) -
                   cv.dens;
          cv.bcons = sqrt_detg * pv.B;
          cv.tracer_ye = cv.dens * pv.ye;
        }
      } else {
        // set to atmo
        cv.bcons(0) = dBx(p.I);
        cv.bcons(1) = dBy(p.I);
        cv.bcons(2) = dBz(p.I);
        pv.B = cv.bcons / sqrt_detg;
        atmo.set(pv, cv, g);
      }
    }

    /* set flag to success */
    con2prim_flag(p.I) = 1;
    // dummy vars
    CCTK_REAL Ex, Ey, Ez, wlor, dumye;

    // Write back pv
    pv.scatter(rho(p.I), eps(p.I), dumye, press(p.I), velx(p.I), vely(p.I),
               velz(p.I), wlor, Ex, Ey, Ez, Bvecx(p.I), Bvecy(p.I), Bvecz(p.I));

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
    if (rep.adjust_cons) {
      cv.scatter(dens(p.I), tau(p.I), dumye, momx(p.I), momy(p.I), momz(p.I),
                 dBx(p.I), dBy(p.I), dBz(p.I));
    }
    // Update saved prims
    saved_rho(p.I) = rho(p.I);
    saved_velx(p.I) = velx(p.I);
    saved_vely(p.I) = vely(p.I);
    saved_velz(p.I) = velz(p.I);
    saved_eps(p.I) = eps(p.I);
  }); // Loop
}

extern "C" void AsterX_Con2Prim_Interpolate_Failed(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim_Interpolate_Failed;
  DECLARE_CCTK_PARAMETERS;

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  const vec<GF3D2<CCTK_REAL>, 6> gf_prims{rho, velx, vely, velz, eps, press};
  const vec<GF3D2<CCTK_REAL>, 5> gf_cons{dens, momx, momy, momz, tau};

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (con2prim_flag(p.I) == 0) {

          const vec<CCTK_REAL, 6> flag_nbs = get_neighbors(con2prim_flag, p);
          const vec<CCTK_REAL, 6> rho_nbs = get_neighbors(rho, p);
          const vec<CCTK_REAL, 6> velx_nbs = get_neighbors(velx, p);
          const vec<CCTK_REAL, 6> vely_nbs = get_neighbors(vely, p);
          const vec<CCTK_REAL, 6> velz_nbs = get_neighbors(velz, p);
          const vec<CCTK_REAL, 6> eps_nbs = get_neighbors(eps, p);
          const vec<CCTK_REAL, 6> saved_rho_nbs = get_neighbors(saved_rho, p);
          const vec<CCTK_REAL, 6> saved_velx_nbs = get_neighbors(saved_velx, p);
          const vec<CCTK_REAL, 6> saved_vely_nbs = get_neighbors(saved_vely, p);
          const vec<CCTK_REAL, 6> saved_velz_nbs = get_neighbors(saved_velz, p);
          const vec<CCTK_REAL, 6> saved_eps_nbs = get_neighbors(saved_eps, p);

          CCTK_REAL sum_nbs =
              sum<6>([&](int i) ARITH_INLINE { return flag_nbs(i); });
          assert(sum_nbs > 0);
          rho(p.I) = calc_avg_neighbors(flag_nbs, rho_nbs, saved_rho_nbs);
          velx(p.I) = calc_avg_neighbors(flag_nbs, velx_nbs, saved_velx_nbs);
          vely(p.I) = calc_avg_neighbors(flag_nbs, vely_nbs, saved_vely_nbs);
          velz(p.I) = calc_avg_neighbors(flag_nbs, velz_nbs, saved_velz_nbs);
          eps(p.I) = calc_avg_neighbors(flag_nbs, eps_nbs, saved_eps_nbs);
          press(p.I) = (gl_gamma - 1) * eps(p.I) * rho(p.I);

          /* reset flag */
          con2prim_flag(p.I) = 1;

          // set to atmos
          /*
          if (rho(p.I) <= rho_abs_min * (1 + atmo_tol))
          {
            const smat<CCTK_REAL, 3> g3_avg([&](int i, int j) ARITH_INLINE
                                            { return calc_avg_v2c(gf_g(i, j),
          p); }); const CCTK_REAL sqrtg = sqrt(calc_det(g3_avg)); const
          vec<CCTK_REAL, 3> Bup{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)}; const
          vec<CCTK_REAL, 3> Blow = calc_contraction(g3_avg, Bup); const
          CCTK_REAL Bsq = calc_contraction(Bup, Blow);

            set_to_atmosphere(rho_abs_min, poly_K, gamma, sqrtg, Bsq, gf_prims,
                              gf_cons, p);
          };
          */
          saved_rho(p.I) = rho(p.I);
          saved_velx(p.I) = velx(p.I);
          saved_vely(p.I) = vely(p.I);
          saved_velz(p.I) = velz(p.I);
          saved_eps(p.I) = eps(p.I);
        }
      });
}

} // namespace AsterX
