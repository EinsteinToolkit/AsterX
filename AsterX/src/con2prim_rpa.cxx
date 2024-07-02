#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include <cmath>

#include "utils.hxx"

#include "reprimand/eos_thermal.h" // The EOS framework
#include "reprimand/eos_idealgas.h"
#include "reprimand/eos_barotropic.h"
#include "reprimand/eos_barotr_poly.h"
#include "reprimand/con2prim_imhd.h" // The con2prim framework

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace EOS_Toolkit;

extern "C" void AsterX_Con2Prim(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim;
  DECLARE_CCTK_PARAMETERS;

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  // Setting up initial data EOS
  const CCTK_REAL n = 1 / (poly_gamma - 1); // Polytropic index
  const CCTK_REAL adiab_ind_id = 1.0 / (poly_gamma - 1);
  const CCTK_REAL rmd_p = pow(poly_k, -n); //Polytropic density scale
  const auto eos_id = make_eos_barotr_poly(adiab_ind_id, rmd_p, rho_max);

  // Setting up evolution EOS
  const CCTK_REAL adiab_ind_evol = 1.0 / (gl_gamma - 1);
  const auto eos = make_eos_idealgas(adiab_ind_evol, eps_max, rho_max);

  // Setting up atmosphere
  const CCTK_REAL rho_atmo_cut = rho_abs_min * (1 + atmo_tol);
  CCTK_REAL eps_atm = eos_id.at_rho(rho_abs_min).eps();
  eps_atm  = eos.range_eps(rho_abs_min, Ye_atmo).limit_to(eps_atm);
  CCTK_REAL p_atm  = eos.at_rho_eps_ye(rho_abs_min, eps_atm, Ye_atmo).press();
  const atmosphere atmo(rho_abs_min, eps_atm, Ye_atmo, p_atm, rho_atmo_cut);

  CCTK_REAL dummy_Ye = 0.5;
  CCTK_REAL dummy_dYe = 0.5;

  // Get a recovery function
  con2prim_mhd cv2pv(eos, 1e-5, 1, 100, 100, atmo, 1e-8, 100);
  //  con2prim_mhd cv2pv(eos, rho_strict, ye_lenient, max_z, max_b, atmo,
  //  c2p_acc,
  //                  max_iter);

  // Loop over the interior of the grid
  cctk_grid.loop_int_device<
      1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                     const PointDesc
                                         &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    /* Get covariant metric */
    const smat<CCTK_REAL, 3> glo(
        [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });

    sm_metric3 g(sm_symt3l(glo(0, 0), glo(0, 1), glo(1, 1), glo(0, 2),
                           glo(1, 2), glo(2, 2)));

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

    cv2pv(pv, cv, g, rep);

    // Handle incorrectable errors
    if (rep.failed()) {
      CCTK_WARN(1, rep.debug_message().c_str());

      if (debug_mode) {
        // need to fix pv to computed values like pv.rho instead of rho(p.I)
        printf(
            "WARNING: "
            "C2P failed. Printing cons and saved prims before set to "
            "atmo: \n"
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
      //set to atmo
      cv.bcons(0) = dBx(p.I);
      cv.bcons(1) = dBy(p.I);
      cv.bcons(2) = dBz(p.I);
      pv.B = cv.bcons / sqrt_detg;
      atmo.set(pv, cv, g);
    }

    /* set flag to success */
    con2prim_flag(p.I) = 1;
    // dummy vars
    CCTK_REAL Ex, Ey, Ez, wlor, dumye;

    // Write back pv
    pv.scatter(rho(p.I), eps(p.I), dumye, press(p.I), velx(p.I), vely(p.I),
               velz(p.I), wlor, Ex, Ey, Ez, Bvecx(p.I), Bvecy(p.I), Bvecz(p.I));

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
