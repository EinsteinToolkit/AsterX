#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include <cmath>

#include "c2p.hxx"
#include "c2p_1DPalenzuela.hxx"
#include "c2p_2DNoble.hxx"

#include <eos_1p.hxx>
#include <eos_polytropic.hxx>

#include <eos.hxx>
#include <eos_idealgas.hxx>

#include "utils.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace EOSX;
using namespace Con2PrimFactory;

enum class eos_t { IdealGas, Hybrid, Tabulated };
enum class c2p_first_t { Noble, Palenzuela };
enum class c2p_second_t { Noble, Palenzuela };

template <typename EOSIDType, typename EOSType>
void AsterX_Con2Prim_typeEoS(CCTK_ARGUMENTS, EOSIDType &eos_cold,
                             EOSType &eos_th) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim;
  DECLARE_CCTK_PARAMETERS;

  c2p_first_t c2p_fir;
  c2p_second_t c2p_sec;

  if (CCTK_EQUALS(c2p_prime, "Noble")) {
    c2p_fir = c2p_first_t::Noble;
  } else if (CCTK_EQUALS(c2p_prime, "Palenzuela")) {
    c2p_fir = c2p_first_t::Palenzuela;
  } else {
    CCTK_ERROR("Unknown value for parameter \"c2p_prime\"");
  }

  if (CCTK_EQUALS(c2p_second, "Noble")) {
    c2p_sec = c2p_second_t::Noble;
  } else if (CCTK_EQUALS(c2p_second, "Palenzuela")) {
    c2p_sec = c2p_second_t::Palenzuela;
  } else {
    CCTK_ERROR("Unknown value for parameter \"c2p_second\"");
  }

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  // Setting up atmosphere
  const CCTK_REAL rho_atmo_cut = rho_abs_min * (1 + atmo_tol);
  const CCTK_REAL gm1 = eos_cold.gm1_from_valid_rmd(rho_abs_min);
  CCTK_REAL eps_atm = eos_cold.sed_from_valid_gm1(gm1);
  eps_atm = std::min(std::max(eos_th.rgeps.min, eps_atm), eos_th.rgeps.max);
  const CCTK_REAL p_atm =
      eos_th.press_from_valid_rho_eps_ye(rho_abs_min, eps_atm, Ye_atmo);
  atmosphere atmo(rho_abs_min, eps_atm, Ye_atmo, p_atm, rho_atmo_cut);

  // Construct Noble c2p object:
  c2p_2DNoble c2p_Noble(eos_th, atmo, max_iter, c2p_tol, rho_strict, vw_lim,
                        B_lim, Ye_lenient);

  // Construct Palenzuela c2p object:
  c2p_1DPalenzuela c2p_Pal(eos_th, atmo, max_iter, c2p_tol, rho_strict, vw_lim,
                           B_lim, Ye_lenient);

  // Loop over the interior of the grid
  cctk_grid.loop_int_device<
      1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                     const PointDesc
                                         &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    /* Get covariant metric */
    const smat<CCTK_REAL, 3> glo(
        [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });

    /* Calculate inverse of 3-metric */
    const CCTK_REAL spatial_detg = calc_det(glo);
    const CCTK_REAL sqrt_detg = sqrt(spatial_detg);

    vec<CCTK_REAL, 3> v_up{saved_velx(p.I), saved_vely(p.I), saved_velz(p.I)};
    const vec<CCTK_REAL, 3> v_low = calc_contraction(glo, v_up);
    CCTK_REAL wlor = calc_wlorentz(v_low, v_up);

    vec<CCTK_REAL, 3> Bup{dBx(p.I) / sqrt_detg, dBy(p.I) / sqrt_detg,
                          dBz(p.I) / sqrt_detg};

    CCTK_REAL dummy_Ye = 0.5;
    CCTK_REAL dummy_dYe = 0.5;
    prim_vars pv;
    prim_vars pv_seeds{saved_rho(p.I), saved_eps(p.I), dummy_Ye, press(p.I),
                       v_up,           wlor,           Bup};
    // Note that cv are densitized, i.e. they all include sqrt_detg
    cons_vars cv{dens(p.I),
                 {momx(p.I), momy(p.I), momz(p.I)},
                 tau(p.I),
                 dummy_dYe,
                 {dBx(p.I), dBy(p.I), dBz(p.I)}};

    if (dens(p.I) <= sqrt_detg * rho_atmo_cut) {
      cv.dBvec(0) = dBx(p.I); // densitized
      cv.dBvec(1) = dBy(p.I);
      cv.dBvec(2) = dBz(p.I);
      pv.Bvec = cv.dBvec / sqrt_detg;
      atmo.set(pv, cv, glo);
      atmo.set(pv_seeds);
    }

    // Construct error report object:
    c2p_report rep_first;
    c2p_report rep_second;

    // Calling the first C2P
    switch (c2p_fir) {
    case c2p_first_t::Noble: {
      c2p_Noble.solve(eos_th, pv, pv_seeds, cv, glo, rep_first);
      break;
    }
    case c2p_first_t::Palenzuela: {
      c2p_Pal.solve(eos_th, pv, pv_seeds, cv, glo, rep_first);
      break;
    }
    default:
      assert(0);
    }

    if (rep_first.failed()) {
      printf("First C2P failed :( \n");
      rep_first.debug_message();
      printf("Calling the back up C2P.. \n");
      // Calling the second C2P
      switch (c2p_sec) {
      case c2p_second_t::Noble: {
        c2p_Noble.solve(eos_th, pv, pv_seeds, cv, glo, rep_second);
        break;
      }
      case c2p_second_t::Palenzuela: {
        c2p_Pal.solve(eos_th, pv, pv_seeds, cv, glo, rep_second);
        break;
      }
      default:
        assert(0);
      }
    }

    if (rep_first.failed() && rep_second.failed()) {
      printf("Second C2P failed too :( :( \n");
      rep_second.debug_message();
      con2prim_flag(p.I) = 0;
    }

    if (rep_first.set_atmo && rep_second.set_atmo) {
      if (debug_mode) {
        printf(
            "WARNING: \n"
            "C2Ps failed. Printing cons and saved prims before set to "
            "atmo: \n"
            "cctk_iteration = %i \n "
            "x, y, z = %26.16e, %26.16e, %26.16e \n "
            "dens = %26.16e \n tau = %26.16e \n momx = %26.16e \n "
            "momy = %26.16e \n momz = %26.16e \n dBx = %26.16e \n "
            "dBy = %26.16e \n dBz = %26.16e \n "
            "saved_rho = %26.16e \n saved_eps = %26.16e \n press= %26.16e \n "
            "saved_velx = %26.16e \n saved_vely = %26.16e \n saved_velz = "
            "%26.16e \n "
            "Bvecx = %26.16e \n Bvecy = %26.16e \n "
            "Bvecz = %26.16e \n "
            "Avec_x = %26.16e \n Avec_y = %26.16e \n Avec_z = %26.16e \n ",
            cctk_iteration, p.x, p.y, p.z, dens(p.I), tau(p.I), momx(p.I),
            momy(p.I), momz(p.I), dBx(p.I), dBy(p.I), dBz(p.I), pv.rho, pv.eps,
            pv.press, pv.vel(0), pv.vel(1), pv.vel(2), pv.Bvec(0), pv.Bvec(1),
            pv.Bvec(2),
            // rho(p.I), eps(p.I), press(p.I), velx(p.I), vely(p.I),
            // velz(p.I), Bvecx(p.I), Bvecy(p.I), Bvecz(p.I),
            Avec_x(p.I), Avec_y(p.I), Avec_z(p.I));
      }

      // set to atmo
      cv.dBvec(0) = dBx(p.I);
      cv.dBvec(1) = dBy(p.I);
      cv.dBvec(2) = dBz(p.I);
      pv.Bvec = cv.dBvec / sqrt_detg;
      atmo.set(pv, cv, glo);

      // assert(0);
    }

    /* set flag to success */
    con2prim_flag(p.I) = 1;

    // dummy vars
    CCTK_REAL Ex, Ey, Ez;

    // Write back pv
    pv.scatter(rho(p.I), eps(p.I), dummy_Ye, press(p.I), velx(p.I), vely(p.I),
               velz(p.I), wlor, Bvecx(p.I), Bvecy(p.I), Bvecz(p.I), Ex, Ey, Ez);

    // Write back cv
    cv.scatter(dens(p.I), momx(p.I), momy(p.I), momz(p.I), tau(p.I), dummy_Ye,
               dBx(p.I), dBy(p.I), dBz(p.I));

    // Update saved prims
    saved_rho(p.I) = rho(p.I);
    saved_velx(p.I) = velx(p.I);
    saved_vely(p.I) = vely(p.I);
    saved_velz(p.I) = velz(p.I);
    saved_eps(p.I) = eps(p.I);
  }); // Loop
}

extern "C" void AsterX_Con2Prim(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Con2Prim;
  DECLARE_CCTK_PARAMETERS;

  // defining EOS objects
  eos_t eostype;
  eos::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
      rgye(ye_min, ye_max);

  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eostype = eos_t::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eostype = eos_t::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated")) {
    eostype = eos_t::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  switch (eostype) {
  case eos_t::IdealGas: {
    CCTK_REAL n = 1 / (poly_gamma - 1); // Polytropic index
    CCTK_REAL rmd_p = pow(poly_k, -n);  // Polytropic density scale

    const eos_polytrope eos_cold(n, rmd_p, rho_max);
    const eos_idealgas eos_th(gl_gamma, particle_mass, rgeps, rgrho, rgye);

    AsterX_Con2Prim_typeEoS(CCTK_PASS_CTOC, eos_cold, eos_th);
    break;
  }
  case eos_t::Hybrid: {
    CCTK_ERROR("Hybrid EOS is not yet supported");
    break;
  }
  case eos_t::Tabulated: {
    CCTK_ERROR("Tabulated EOS is not yet supported");
    break;
  }
  default:
    assert(0);
  }
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
