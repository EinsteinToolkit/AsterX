#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "c2p.hxx"
#include "c2p_1DEntropy.hxx"
#include "c2p_1DPalenzuela.hxx"
#include "c2p_1DRePrimAnd.hxx"
#include "c2p_2DNoble.hxx"

#include "c2p_utils.hxx"

#include "setup_eos.hxx"

#include <array>

namespace Con2PrimFactory {

using namespace Arith;
using namespace EOSX;
using namespace AsterUtils;

namespace {

void check_pal(const char *name, CCTK_REAL actual, CCTK_REAL expected) {
  if (!std::isfinite(actual) || !std::isfinite(expected) ||
      fabs(actual - expected) > 1.0e-8 * fmax(fabs(expected), 1.0e-6))
    CCTK_VERROR("Con2PrimFactory test %s: actual=%.16e expected=%.16e",
               name, actual, expected);
}

void check_cons(const cons_vars &cv, const cons_vars &expected) {
  check_pal("dens", cv.dens, expected.dens);
  check_pal("tau", cv.tau, expected.tau);
  check_pal("DYe", cv.DYe, expected.DYe);
  check_pal("DEnt", cv.DEnt, expected.DEnt);
  for (int d = 0; d < 3; ++d) {
    check_pal("mom", cv.mom(d), expected.mom(d));
    check_pal("dBvec", cv.dBvec(d), expected.dBvec(d));
  }
}

template <typename EOSType>
void test_pal(const EOSType &eos, bool use_temp) {
  const CCTK_REAL rho_atmo = eos.rgrho.min;
  const CCTK_REAL Ye_atmo = eos.rgye.max;
  const CCTK_REAL temp_atmo = eos.rgtemp.min;
  CCTK_REAL eps_atmo =
      eos.eps_from_rho_temp_ye(rho_atmo, temp_atmo, Ye_atmo);
  const CCTK_REAL press_atmo =
      eos.press_from_rho_temp_ye(rho_atmo, temp_atmo, Ye_atmo);
  const CCTK_REAL entropy_atmo =
      eos.kappa_from_rho_eps_ye(rho_atmo, eps_atmo, Ye_atmo);
  atmosphere atmo(rho_atmo, eps_atmo, Ye_atmo, press_atmo, temp_atmo,
                  entropy_atmo, rho_atmo * 1.001);
  c2p_1DPalenzuela c2p_Pal(&eos, atmo, 200, 1.0e-12, -1.0, 10.0, 100.0,
                           1.0e20, 1.0e20, 1.0e20, 1.0e20, 1.0e20,
                           true, false, use_temp, false, false, 1.0);
  // Non-unit determinant also checks the densitized conservative update.
  const smat<CCTK_REAL, 3> g{1.2, 0.0, 0.0, 1.1, 0.0, 0.9};
  const vec<CCTK_REAL, 3> beta{0.0, 0.0, 0.0};

  for (CCTK_REAL rho : {3.0e-4, 3.0e-3})
    for (CCTK_REAL Ye : {0.15, 0.45})
      for (CCTK_REAL temp : {0.0011, 0.02, 0.07})
        for (CCTK_REAL v : {0.0, 0.2})
          for (CCTK_REAL B : {0.0, 0.001}) {
            CCTK_REAL eps = eos.eps_from_rho_temp_ye(rho, temp, Ye);
            const CCTK_REAL press = eos.press_from_rho_temp_ye(rho, temp, Ye);
            const CCTK_REAL entropy = eos.kappa_from_rho_eps_ye(rho, eps, Ye);
            const vec<CCTK_REAL, 3> vel{v, 0.0, 0.0};
            const CCTK_REAL wlor = calc_wlorentz(calc_contraction(g, vel), vel);
            const prim_vars pv_in{rho, eps, Ye, press, temp, entropy,
                                  vel, wlor, {B, 0.2 * B, 0.0}};
            cons_vars cv_in;
            cv_in.from_prim(pv_in, g);
            for (bool reject : {false, true}) {
              prim_vars pv;
              cons_vars cv = cv_in;
              c2p_report rep;
              c2p_Pal.solve(&eos, pv, cv, 1.0, beta, g, rep, reject);
              if (rep.failed() || rep.set_atmo || rep.adjust_cons)
                CCTK_ERROR("Palenzuela test: valid interior state was rejected "
                           "or adjusted");
              check_pal("rho", pv.rho, rho);
              check_pal("eps", pv.eps, eps);
              check_pal("Ye", pv.Ye, Ye);
              check_pal("temperature", pv.temperature, temp);
              check_pal("press", pv.press, press);
              check_pal("entropy", pv.entropy, entropy);
              check_pal("wlor", pv.w_lor, wlor);
              for (int d = 0; d < 3; ++d) {
                check_pal("vel", pv.vel(d), vel(d));
                check_pal("Bvec", pv.Bvec(d), pv_in.Bvec(d));
              }
              check_cons(cv, cv_in);
            }
          }

  if (!use_temp)
    return;

  // Stationary, unmagnetized input gives eps_raw = tau / D exactly.
  // Test both bounds and the old ideal-gas zero-energy fallback policy.
  const CCTK_REAL rho = 3.0e-4, Ye = 0.3;
  const auto rgeps = eos.range_eps_from_rho_ye(rho, Ye);
  const CCTK_REAL dens = sqrt(calc_det(g)) * rho;
  for (CCTK_REAL eps_raw : {rgeps.min - 1.0e-4, 0.0, rgeps.max + 1.0e-4})
    for (bool reject : {false, true}) {
      const cons_vars cv_in{dens, {0.0, 0.0, 0.0}, dens * eps_raw,
                            dens * Ye, 0.0, {0.0, 0.0, 0.0}};
      cons_vars cv = cv_in;
      prim_vars pv;
      c2p_report rep;
      c2p_Pal.solve(&eos, pv, cv, 1.0, beta, g, rep, reject);
      const bool expect_fail =
          reject && (rgeps.min < 0.0 ? eps_raw < rgeps.min : eps_raw <= 0.0);
      if (expect_fail) {
        if (rep.status != c2p_report::RANGE_EPS)
          CCTK_ERROR("Palenzuela test: raw energy did not trigger fallback");
        check_cons(cv, cv_in);
        continue;
      }
      const bool clipped = eps_raw < rgeps.min || eps_raw > rgeps.max;
      if (rep.failed() || rep.set_atmo || rep.adjust_cons != clipped)
        CCTK_ERROR("Palenzuela test: incorrect energy-clipping report");
      CCTK_REAL eps = std::min(std::max(eps_raw, rgeps.min), rgeps.max);
      check_pal("bounded eps", pv.eps, eps);
      check_pal("bounded temperature", pv.temperature,
                eos.temp_from_rho_eps_ye(rho, eps, Ye));
      check_pal("bounded pressure", pv.press,
                eos.press_from_rho_eps_ye(rho, eps, Ye));
      check_pal("bounded entropy", pv.entropy,
                eos.kappa_from_rho_eps_ye(rho, eps, Ye));
      cons_vars expected;
      expected.from_prim(pv, g);
      check_cons(cv, expected);
    }
}

template <typename EOSType>
void test_rpa(const EOSType &eos, bool use_temp) {
  const CCTK_REAL rho_atmo = eos.rgrho.min;
  const CCTK_REAL Ye_atmo = eos.rgye.max;
  const CCTK_REAL temp_atmo = eos.rgtemp.min;
  CCTK_REAL eps_atmo =
      eos.eps_from_rho_temp_ye(rho_atmo, temp_atmo, Ye_atmo);
  const CCTK_REAL press_atmo =
      eos.press_from_rho_temp_ye(rho_atmo, temp_atmo, Ye_atmo);
  const CCTK_REAL entropy_atmo =
      eos.kappa_from_rho_eps_ye(rho_atmo, eps_atmo, Ye_atmo);
  atmosphere atmo(rho_atmo, eps_atmo, Ye_atmo, press_atmo, temp_atmo,
                  entropy_atmo, rho_atmo * 1.001);
  c2p_1DRePrimAnd c2p_RPA(&eos, atmo, 200, 1.0e-12, -1.0, 10.0, 100.0,
                          1.0e20, 1.0e20, 1.0e20, 1.0e20, 1.0e20,
                          true, false, use_temp, false, false, 1.0);
  const smat<CCTK_REAL, 3> g{1.2, 0.0, 0.0, 1.1, 0.0, 0.9};
  const smat<CCTK_REAL, 3> gup = calc_inv(g, calc_det(g));
  const CCTK_REAL sqrt_detg = sqrt(calc_det(g));
  const vec<CCTK_REAL, 3> beta{0.0, 0.0, 0.0};

  // Reuse the scalar/conservative checks from the Palenzuela tests.
  // The cold shifted-table states include h < 1 and therefore mu > 1 at rest.
  for (CCTK_REAL rho : {3.0e-4, 3.0e-3})
    for (CCTK_REAL Ye : {0.15, 0.45})
      for (CCTK_REAL temp : {0.0011, 0.02, 0.07})
        for (CCTK_REAL v : {0.0, 0.2, 0.75})
          for (CCTK_REAL B : {0.0, 0.001}) {
            CCTK_REAL eps = eos.eps_from_rho_temp_ye(rho, temp, Ye);
            const CCTK_REAL press = eos.press_from_rho_temp_ye(rho, temp, Ye);
            const CCTK_REAL entropy = eos.kappa_from_rho_eps_ye(rho, eps, Ye);
            const vec<CCTK_REAL, 3> vel{v, 0.0, 0.0};
            const CCTK_REAL wlor = calc_wlorentz(vel, calc_contraction(g, vel));
            const prim_vars pv_in{rho, eps, Ye, press, temp, entropy,
                                  vel, wlor, {B, 0.2 * B, 0.0}};
            cons_vars cv_in;
            cv_in.from_prim(pv_in, g);

            const CCTK_REAL d = rho * wlor;
            const vec<CCTK_REAL, 3> r = cv_in.mom / cv_in.dens;
            const vec<CCTK_REAL, 3> b = pv_in.Bvec / sqrt(d);
            const CCTK_REAL rb = calc_contraction(r, b);
            typename RePrimAnd::froot<EOSType>::cache cache{};
            RePrimAnd::froot<EOSType> f(
                &eos, Ye, d, cv_in.tau / cv_in.dens,
                calc_contraction(calc_contraction(gup, r), r), rb * rb,
                calc_contraction(calc_contraction(g, b), b), cache);
            const CCTK_REAL h = 1.0 + eps + press / rho;
            if (!(f.h0 > 0.0 && f.h0 <= h))
              CCTK_ERROR("RePrimAnd test: invalid enthalpy lower bound");
            ROOTSTAT status{};
            const auto bracket = f.initial_bracket(status);
            const CCTK_REAL mu = 1.0 / (h * wlor);
            if (status != ROOTSTAT::SUCCESS ||
                !(bracket.min() <= mu && mu <= bracket.max()))
              CCTK_ERROR("RePrimAnd test: physical root is outside the bracket");
            check_pal("RPA master function at physical root", f(mu), 0.0);

            for (bool reject : {false, true}) {
              prim_vars pv;
              cons_vars cv = cv_in;
              c2p_report rep;
              c2p_RPA.solve(&eos, pv, cv, 1.0, beta, g, rep, reject);
              if (rep.failed() || rep.set_atmo || rep.adjust_cons)
                CCTK_ERROR("RePrimAnd test: valid interior state was rejected "
                           "or adjusted");
              check_pal("RPA rho", pv.rho, rho);
              check_pal("RPA eps", pv.eps, eps);
              check_pal("RPA Ye", pv.Ye, Ye);
              check_pal("RPA temperature", pv.temperature, temp);
              check_pal("RPA press", pv.press, press);
              check_pal("RPA entropy", pv.entropy, entropy);
              check_pal("RPA wlor", pv.w_lor, wlor);
              for (int d = 0; d < 3; ++d) {
                check_pal("RPA vel", pv.vel(d), vel(d));
                check_pal("RPA Bvec", pv.Bvec(d), pv_in.Bvec(d));
              }
              check_cons(cv, cv_in);
            }
          }

  const CCTK_REAL rho = 3.0e-4, Ye = 0.3;
  const auto rgeps = eos.range_eps_from_rho_ye(rho, Ye);
  const CCTK_REAL dens = sqrt_detg * rho;
  // RePrimAnd already clips local energy for both choices of use_temp.
  for (CCTK_REAL eps_raw : {rgeps.min - 1.0e-4, 0.0, rgeps.max + 1.0e-4})
    for (bool reject : {false, true}) {
      const cons_vars cv_in{dens, {0.0, 0.0, 0.0}, dens * eps_raw,
                            dens * Ye, 0.0, {0.0, 0.0, 0.0}};
      cons_vars cv = cv_in;
      prim_vars pv;
      c2p_report rep;
      c2p_RPA.solve(&eos, pv, cv, 1.0, beta, g, rep, reject);
      const bool expect_fail =
          reject && (rgeps.min < 0.0 ? eps_raw < rgeps.min : eps_raw <= 0.0);
      if (expect_fail) {
        if (rep.status != c2p_report::RANGE_EPS)
          CCTK_ERROR("RePrimAnd test: raw energy did not trigger fallback");
        check_cons(cv, cv_in);
        continue;
      }
      const bool clipped = eps_raw < rgeps.min || eps_raw > rgeps.max;
      if (rep.failed() || rep.set_atmo || rep.adjust_cons != clipped)
        CCTK_ERROR("RePrimAnd test: incorrect energy-clipping report");
      CCTK_REAL eps = std::min(std::max(eps_raw, rgeps.min), rgeps.max);
      check_pal("RPA bounded eps", pv.eps, eps);
      check_pal("RPA bounded temperature", pv.temperature,
                eos.temp_from_rho_eps_ye(rho, eps, Ye));
      check_pal("RPA bounded pressure", pv.press,
                eos.press_from_rho_eps_ye(rho, eps, Ye));
      check_pal("RPA bounded entropy", pv.entropy,
                eos.kappa_from_rho_eps_ye(rho, eps, Ye));
      cons_vars expected;
      expected.from_prim(pv, g);
      check_cons(cv, expected);
    }

  for (CCTK_REAL tau : {std::numeric_limits<CCTK_REAL>::quiet_NaN(),
                        std::numeric_limits<CCTK_REAL>::infinity(),
                        -std::numeric_limits<CCTK_REAL>::infinity()}) {
    cons_vars cv{dens, {0.0, 0.0, 0.0}, tau,
                  dens * Ye, 0.0, {0.0, 0.0, 0.0}};
    prim_vars pv;
    c2p_report rep;
    c2p_RPA.solve(&eos, pv, cv, 1.0, beta, g, rep);
    if (!rep.failed() ||
        !(std::isnan(tau) ? std::isnan(cv.tau) : cv.tau == tau))
      CCTK_ERROR("RePrimAnd test: non-finite energy was accepted or overwritten");
    cv.tau = 0.0;
    const cons_vars expected{dens, {0.0, 0.0, 0.0}, 0.0,
                             dens * Ye, 0.0, {0.0, 0.0, 0.0}};
    check_cons(cv, expected);
  }

  // Invalid EOS bounds must fail before constructing 1/h0.
  for (CCTK_REAL eps_min :
       {-1.0, -2.0, std::numeric_limits<CCTK_REAL>::quiet_NaN()}) {
    auto bad = eos;
    bad.rgeps.min = eps_min;
    const cons_vars cv_in{dens, {0.0, 0.0, 0.0}, 0.0,
                          dens * Ye, 0.0, {0.0, 0.0, 0.0}};
    cons_vars cv = cv_in;
    prim_vars pv;
    c2p_report rep;
    c2p_RPA.solve(&bad, pv, cv, 1.0, beta, g, rep);
    if (rep.status != c2p_report::RANGE_EPS)
      CCTK_ERROR("RePrimAnd test: unsupported EOS energy range was accepted");
    check_cons(cv, cv_in);
  }
}

template <typename EOSType>
void test_cons(const EOSType &eos) {
  // This limiter uses no solver or atmosphere state.
  c2p c2p_test{};
  const smat<CCTK_REAL, 3> g{1.2, 0.0, 0.0, 1.1, 0.0, 0.9};
  const CCTK_REAL sqrt_detg = sqrt(calc_det(g));
  const CCTK_REAL tauFluid_atmo = 1.0e-15;

  // Physical states, including negative table eps, must remain unchanged.
  for (CCTK_REAL rho : {3.0e-4, 3.0e-3})
    for (CCTK_REAL Ye : {0.15, 0.45})
      for (CCTK_REAL temp : {0.0011, 0.02, 0.07})
        for (CCTK_REAL v : {0.0, 0.2, 0.75})
          for (CCTK_REAL B : {0.0, 0.001}) {
            CCTK_REAL eps = eos.eps_from_rho_temp_ye(rho, temp, Ye);
            const CCTK_REAL press = eos.press_from_rho_temp_ye(rho, temp, Ye);
            const CCTK_REAL entropy = eos.kappa_from_rho_eps_ye(rho, eps, Ye);
            const vec<CCTK_REAL, 3> vel{v, 0.0, 0.0};
            const CCTK_REAL wlor = calc_wlorentz(vel, calc_contraction(g, vel));
            const prim_vars pv{rho, eps, Ye, press, temp, entropy,
                               vel, wlor, {B, 0.2 * B, 0.0}};
            cons_vars cv_in;
            cv_in.from_prim(pv, g);
            cons_vars cv = cv_in;
            c2p_test.cons_floors_and_ceilings(&eos, cv, g, tauFluid_atmo);
            if (cv.tau != cv_in.tau)
              CCTK_ERROR("Conservative test: valid energy was modified");
            check_cons(cv, cv_in);
          }

  // Exercise the trigger and repair separately, including density/Ye bounds.
  // Equality is not repaired; no margin is added to an already valid bound.
  for (CCTK_REAL rho : {0.5 * eos.rgrho.min, 3.0e-4, 2.0 * eos.rgrho.max})
    for (CCTK_REAL Ye : {0.0, 0.3, 1.0})
      for (CCTK_REAL B : {0.0, 0.001}) {
        const CCTK_REAL dens = sqrt_detg * rho;
        const vec<CCTK_REAL, 3> dBvec =
            sqrt_detg * vec<CCTK_REAL, 3>{B, 0.2 * B, 0.0};
        const CCTK_REAL tau_mag =
            0.5 * calc_contraction(calc_contraction(g, dBvec), dBvec) /
            sqrt_detg;
        const CCTK_REAL tau_lim =
            tau_mag + dens * fmin(0.0, eos.rgeps.min);
        const CCTK_REAL rhoL = fmin(fmax(rho, eos.rgrho.min), eos.rgrho.max);
        const CCTK_REAL YeL = fmin(fmax(Ye, eos.rgye.min), eos.rgye.max);
        const auto rgeps = eos.range_eps_from_rho_ye(rhoL, YeL);
        for (CCTK_REAL offset : {-1.0e-4, 0.0, 1.0e-4}) {
          const cons_vars cv_in{dens, {0.0, 0.0, 0.0},
                                tau_lim + dens * offset, dens * Ye, 0.0,
                                dBvec};
          cons_vars cv = cv_in, expected = cv_in;
          if (offset < 0.0)
            expected.tau =
                tau_mag + dens * rgeps.min + sqrt_detg * tauFluid_atmo;
          c2p_test.cons_floors_and_ceilings(&eos, cv, g, tauFluid_atmo);
          if (offset >= 0.0 && cv.tau != cv_in.tau)
            CCTK_ERROR("Conservative test: energy at or above bound changed");
          check_cons(cv, expected);
        }
      }

  // Invalid density must not be used in DYe / D or an EOS query.
  for (CCTK_REAL dens : {0.0, -1.0}) {
    const cons_vars cv_in{dens, {0.0, 0.0, 0.0}, -1.0, 0.0, 0.0,
                          {0.0, 0.0, 0.0}};
    cons_vars cv = cv_in;
    c2p_test.cons_floors_and_ceilings(&eos, cv, g, tauFluid_atmo);
    check_cons(cv, cv_in);
  }

  // Do not turn a non-finite energy into an apparently valid state.
  for (CCTK_REAL tau : {std::numeric_limits<CCTK_REAL>::quiet_NaN(),
                        std::numeric_limits<CCTK_REAL>::infinity(),
                        -std::numeric_limits<CCTK_REAL>::infinity()}) {
    cons_vars cv{1.0, {0.0, 0.0, 0.0}, tau, 0.3, 0.0, {0.0, 0.0, 0.0}};
    c2p_test.cons_floors_and_ceilings(&eos, cv, g, tauFluid_atmo);
    if (!(std::isnan(tau) ? std::isnan(cv.tau) : cv.tau == tau))
      CCTK_ERROR("Conservative test: non-finite energy was overwritten");
    cv.tau = 0.0;
    const cons_vars expected{1.0, {0.0, 0.0, 0.0}, 0.0, 0.3, 0.0,
                             {0.0, 0.0, 0.0}};
    check_cons(cv, expected);
  }

  // The existing momentum cap must still enforce |S| <= D + tau.
  cons_vars cv{1.0, {10.0, 0.0, 0.0}, 1.0, 0.3, 0.0, {0.0, 0.0, 0.0}};
  c2p_test.cons_floors_and_ceilings(&eos, cv, g, tauFluid_atmo);
  const cons_vars expected{1.0, {2.0 * sqrt(g(0, 0)), 0.0, 0.0}, 1.0,
                           0.3, 0.0, {0.0, 0.0, 0.0}};
  check_cons(cv, expected);
}

void test_atmo_reset(const atmosphere &atmo) {
  // A non-diagonal metric also checks the magnetic energy and densitization.
  const smat<CCTK_REAL, 3> g{1.2, 0.1, 0.0, 1.1, 0.05, 0.9};
  const CCTK_REAL sqrt_detg = sqrt(calc_det(g));

  for (CCTK_REAL B : {0.0, 0.001}) {
    const vec<CCTK_REAL, 3> Bup{B, 0.2 * B, -0.1 * B};
    prim_vars pv;
    pv.set_to_nan();
    pv.Bvec = Bup;

    // The reset must overwrite every fluid field, and be idempotent.
    for (int repeat = 0; repeat < 2; ++repeat) {
      atmo.set(pv);
      if (pv.rho != atmo.rho_atmo || pv.eps != atmo.eps_atmo ||
          pv.Ye != atmo.ye_atmo || pv.press != atmo.press_atmo ||
          pv.temperature != atmo.temp_atmo ||
          pv.entropy != atmo.entropy_atmo || pv.w_lor != 1.0)
        CCTK_ERROR("Atmosphere reset did not copy the complete thermal state");
      for (int d = 0; d < 3; ++d)
        if (pv.vel(d) != 0.0 || pv.E(d) != 0.0 || pv.Bvec(d) != Bup(d))
          CCTK_ERROR("Atmosphere reset changed B or retained velocity/E");
    }

    const CCTK_REAL dens = sqrt_detg * atmo.rho_atmo;
    const CCTK_REAL Bsq = calc_contraction(Bup, calc_contraction(g, Bup));
    const cons_vars expected{
        dens, {0.0, 0.0, 0.0},
        dens * atmo.eps_atmo + 0.5 * sqrt_detg * Bsq,
        dens * atmo.ye_atmo, dens * atmo.entropy_atmo, sqrt_detg * Bup};

    // Rebuilding conservatives must agree with the stationary atmosphere.
    cons_vars cv;
    cv.from_prim(pv, g);
    check_cons(cv, expected);

    // The joint reset must replace dirty conservatives and preserve dB.
    pv.set_to_nan();
    pv.Bvec = Bup;
    cv.set_to_nan();
    cv.dBvec = sqrt_detg * Bup;
    atmo.set(pv, cv, g);
    check_cons(cv, expected);
    if (pv.rho != atmo.rho_atmo || pv.eps != atmo.eps_atmo ||
        pv.Ye != atmo.ye_atmo || pv.press != atmo.press_atmo ||
        pv.temperature != atmo.temp_atmo ||
        pv.entropy != atmo.entropy_atmo || pv.w_lor != 1.0)
      CCTK_ERROR("Joint atmosphere reset left an inconsistent primitive state");
    for (int d = 0; d < 3; ++d)
      if (pv.vel(d) != 0.0 || pv.E(d) != 0.0 || pv.Bvec(d) != Bup(d))
        CCTK_ERROR("Joint atmosphere reset changed B or retained velocity/E");
  }
}

template <typename EOSType>
void check_atmo(const atmosphere &atmo, const EOSType &eos,
                CCTK_REAL rho, CCTK_REAL temp, CCTK_REAL Ye,
                CCTK_REAL atmo_tol) {
  check_pal("atmo rho", atmo.rho_atmo, rho);
  check_pal("atmo temperature", atmo.temp_atmo, temp);
  check_pal("atmo Ye", atmo.ye_atmo, Ye);
  check_pal("atmo eps", atmo.eps_atmo,
            eos.eps_from_rho_temp_ye(rho, temp, Ye));
  check_pal("atmo pressure", atmo.press_atmo,
            eos.press_from_rho_temp_ye(rho, temp, Ye));
  check_pal("atmo kappa", atmo.entropy_atmo,
            eos.kappa_from_rho_temp_ye(rho, temp, Ye));
  check_pal("atmo cutoff", atmo.rho_cut, rho * (1 + atmo_tol));
  test_atmo_reset(atmo);
}

template <typename EOSType>
void test_atmo(const EOSType &eos) {
  // Thermal atmosphere must not access the cold EOS.
  const eos_1p_polytropic *eos_1p = nullptr;
  const CCTK_REAL r_atmo = 10.0, atmo_tol = 0.001;
  const CCTK_REAL rho_abs_min = 1.0e-3, t_atmo = 0.02, Ye_atmo = 0.3;

  // Constant, graded, inner/boundary/outer states, and distinct neighbors.
  for (CCTK_REAL r : {0.0, 5.0, 10.0, 15.0, 40.0, 100.0})
    for (CCTK_REAL nr : {0.0, 6.0})
      for (CCTK_REAL nt : {0.0, 2.0}) {
        const CCTK_REAL f = r > r_atmo ? r_atmo / r : 1.0;
        const CCTK_REAL rho =
            std::clamp(rho_abs_min * pow(f, nr), eos.rgrho.min, eos.rgrho.max);
        const CCTK_REAL temp =
            std::clamp(t_atmo * pow(f, nt), eos.rgtemp.min, eos.rgtemp.max);
        // Pressure and its exponent are inactive in temperature-primary mode.
        for (CCTK_REAL p_atmo : {0.0, 1.0})
          for (CCTK_REAL np : {0.0, 12.0}) {
            const auto atmo = make_atmo(
                eos_1p, &eos, r, rho_abs_min, p_atmo, t_atmo, Ye_atmo,
                r_atmo, nr, np, nt, atmo_tol, true, false, false);
            check_atmo(atmo, eos, rho, temp, Ye_atmo, atmo_tol);
          }
      }

  // Inputs at/beyond the EOS bounds must produce a complete bounded state.
  for (CCTK_REAL rho : {0.0, 0.5 * eos.rgrho.min, 2.0 * eos.rgrho.max})
    for (CCTK_REAL temp : {0.5 * eos.rgtemp.min, 2.0 * eos.rgtemp.max})
      for (CCTK_REAL Ye : {eos.rgye.min - 0.1, eos.rgye.max + 0.1})
        for (CCTK_REAL tol : {0.0, 0.001, 0.1}) {
          const auto atmo =
              make_atmo(eos_1p, &eos, 0.0, rho, 0.0, temp, Ye, r_atmo, 6.0,
                        12.0, 2.0, tol, true, false, false);
          check_atmo(atmo, eos,
                      std::clamp(rho, eos.rgrho.min, eos.rgrho.max),
                      std::clamp(temp, eos.rgtemp.min, eos.rgtemp.max),
                      std::clamp(Ye, eos.rgye.min, eos.rgye.max), tol);
        }
}

void test_atmo_beq(const eos_3p_tabulated3d &eos) {
  const eos_1p_polytropic *eos_1p = nullptr;
  const CCTK_REAL r_atmo = 10.0, rho_abs_min = 1.0e-3;
  const CCTK_REAL t_atmo = 0.02, atmo_tol = 0.001;

  // The synthetic chemical potentials set a rho- and T-dependent root.
  for (CCTK_REAL r : {0.0, 10.0, 20.0, 100.0}) {
    const CCTK_REAL f = r > r_atmo ? r_atmo / r : 1.0;
    const CCTK_REAL rho =
        std::clamp(rho_abs_min * pow(f, 6.0), eos.rgrho.min, eos.rgrho.max);
    const CCTK_REAL temp =
        std::clamp(t_atmo * pow(f, 2.0), eos.rgtemp.min, eos.rgtemp.max);
    const CCTK_REAL Ye =
        0.3 + 0.01 * (log(rho) - log(1.0e-4)) +
        0.02 * (log(temp) - log(1.0e-2));
    const auto atmo =
        make_atmo(eos_1p, &eos, r, rho_abs_min, 0.0, t_atmo, 0.49,
                  r_atmo, 6.0, 0.0, 2.0, atmo_tol, true, false, true);
    check_atmo(atmo, eos, rho, temp, Ye, atmo_tol);
  }
}

void test_atmo_ideal(const eos_3p_idealgas &eos_in) {
  // P_cold = 100 rho^2, eps_cold = 100 rho.
  eos_1p_polytropic eos_1p;
  eos_1p.init(2.0, 100.0, eos_in.rgrho.max);
  const CCTK_REAL r_atmo = 10.0, atmo_tol = 0.001;
  const CCTK_REAL rho_abs_min = 1.0e-3, Ye_atmo = 0.3;

  // Check particle-mass conversion, also when cold/evolution gamma differ.
  for (CCTK_REAL umass : {1.0, 2.0}) {
    eos_3p_idealgas eos;
    auto rgeps = eos_in.rgeps;
    eos.init(eos_in.gamma, umass, rgeps, eos_in.rgrho, eos_in.rgye);
    for (CCTK_REAL r : {0.0, 10.0, 20.0, 40.0, 100.0})
      for (CCTK_REAL p_atmo : {0.0, 1.0e-6, 1.0})
        for (CCTK_REAL np : {0.0, 12.0})
          for (bool thermal : {false, true})
            for (bool use_press : {false, true}) {
              if (thermal && !use_press)
                continue; // Temperature-primary mode is tested above.
              const CCTK_REAL f = r > r_atmo ? r_atmo / r : 1.0;
              const CCTK_REAL rho =
                  std::clamp(rho_abs_min * pow(f, 6.0),
                             eos.rgrho.min, eos.rgrho.max);
              CCTK_REAL eps = thermal
                  ? p_atmo * pow(f, np) / ((eos.gamma - 1.0) * rho)
                  : 100.0 * rho;
              eps = std::clamp(eps, eos.rgeps.min, eos.rgeps.max);
              const CCTK_REAL temp = (eos.gamma - 1.0) * umass * eps;

              // Temperature and its exponent are inactive in these modes.
              // Cold matching also ignores p_atmo, np and use_press.
              const auto atmo = make_atmo(
                  &eos_1p, &eos, r, rho_abs_min, p_atmo,
                  2.0 * eos.rgtemp.max, Ye_atmo, r_atmo, 6.0, np, 7.0,
                  atmo_tol, thermal, use_press, false);
              check_atmo(atmo, eos, rho, temp, Ye_atmo, atmo_tol);
            }
  }
}

void test_pal_energy() {
  // Host-local synthetic table: exercise the actual table inverse without
  // loading a production EOS. No host pointers are captured in GPU kernels.
  std::array<CCTK_REAL, 3> lr{log(1.0e-6), log(1.0e-4), log(1.0e-2)};
  std::array<CCTK_REAL, 3> lt{log(1.0e-3), log(1.0e-2), log(1.0e-1)};
  std::array<CCTK_REAL, 2> ye{0.1, 0.5};
  std::array<CCTK_REAL, 18 * NTABLES> data{};
  for (int k = 0; k < 2; ++k)
    for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 3; ++i) {
        const int offset = NTABLES * (i + 3 * (j + 3 * k));
        data[offset + eos_3p_tabulated3d::PRESS] = lr[i] + lt[j];
        // Cold dense states can have eps below eps_atmo at the same T.
        data[offset + eos_3p_tabulated3d::EPS] =
            lt[j] - 0.1 * lr[i] + 0.2 * ye[k];
        data[offset + eos_3p_tabulated3d::S] = lt[j] - 0.1 * lr[i] + ye[k];
        const CCTK_REAL Ye_beq =
            0.3 + 0.01 * (lr[i] - lr[1]) + 0.02 * (lt[j] - lt[1]);
        data[offset + eos_3p_tabulated3d::MU_E] = ye[k] - Ye_beq;
        data[offset + eos_3p_tabulated3d::MU_P] = 0.0;
        data[offset + eos_3p_tabulated3d::MU_N] = 0.0;
      }
  linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> interp(
      data.data(), {3, 3, 2}, lr.data(), lt.data(), ye.data());
  for (CCTK_REAL shift : {0.0, 0.05}) {
    eos_3p_tabulated3d eos;
    eos.gamma = 1.4; // Unused by Palenzuela's root equation.
    eos.interptable = &interp;
    eos.energy_shift = &shift;
    eos.rgrho = {exp(lr.front()), exp(lr.back())};
    eos.rgtemp = {exp(lt.front()), exp(lt.back())};
    eos.rgye = {ye.front(), ye.back()};
    eos.rgeps = eos.compute_eps_range_full_table();
    test_atmo(eos);
    test_atmo_beq(eos);
    test_cons(eos);
    test_pal(eos, true);
    test_rpa(eos, true);
  }
  for (CCTK_REAL gamma : {1.4, 2.0}) {
    eos_3p_idealgas eos;
    eos_3p::range er{0.0, 1.0}, rr{1.0e-6, 1.0e-2}, yr{0.1, 0.5};
    eos.init(gamma, 1.0, er, rr, yr);
    test_atmo(eos);
    test_atmo_ideal(eos);
    test_cons(eos);
    test_pal(eos, false);
    test_pal(eos, true);
    test_rpa(eos, false);
    test_rpa(eos, true);

    // A positive ideal-gas eps_min must not become a new rejection policy.
    er.min = 0.001;
    eos.init(gamma, 1.0, er, rr, yr);
    test_atmo(eos);
    test_atmo_ideal(eos);
    test_cons(eos);
    test_rpa(eos, false);
    test_rpa(eos, true);
  }
  CCTK_INFO("Atmosphere, conservative and C2P energy-bound tests passed");
}

} // namespace

extern "C" void Con2PrimFactory_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  test_pal_energy();

  // The existing multi-solver example below requires an active ideal-gas EOS.
  auto eos_3p_ig = global_eos_3p_ig;
  if (!eos_3p_ig)
    return;

  // Set atmo values
  const CCTK_REAL rho_atmo = 1e-10;
  CCTK_REAL eps_atmo = 1e-8;
  const CCTK_REAL Ye_atmo = 0.5;
  const CCTK_REAL press_atmo =
      eos_3p_ig->press_from_rho_eps_ye(rho_atmo, eps_atmo, Ye_atmo);
  const CCTK_REAL temp_atmo =
      eos_3p_ig->temp_from_rho_eps_ye(rho_atmo, eps_atmo, Ye_atmo);
  const CCTK_REAL entropy_atmo =
      eos_3p_ig->kappa_from_rho_eps_ye(rho_atmo, eps_atmo, Ye_atmo);

  // Setting up atmosphere
  const CCTK_REAL rho_atmo_cut = rho_atmo * (1 + 1.0e-3);
  atmosphere atmo(rho_atmo, eps_atmo, Ye_atmo, press_atmo, temp_atmo,
                  entropy_atmo, rho_atmo_cut);

  // Metric

  const CCTK_REAL alp = 1.0;

  const vec<CCTK_REAL, 3> beta{0.0, 0.0, 0.0};

  const smat<CCTK_REAL, 3> g{1.0, 0.0, 0.0,
                             1.0, 0.0, 1.0}; // xx, xy, xz, yy, yz, zz

  // Set BH limiters
  const CCTK_REAL alp_thresh = -1.;
  const CCTK_REAL rho_BH = 1e20;
  const CCTK_REAL eps_BH = 1e20;
  const CCTK_REAL vwlim_BH = 1e20;

  // Mag. limits
  const CCTK_REAL sigma_max = 100.;
  const CCTK_REAL inv_beta_max = 100.;

  // Con2Prim objects
  // (eos_3p_ig, atmo, max_iter, c2p_tol
  //  alp_thresh,
  //  vw_lim, B_lim, rho_BH, eps_BH, vwlim_BH,
  //  Ye_lenient, use_z, use_temperature)
  c2p_2DNoble c2p_Noble(eos_3p_ig, atmo, 100, 1e-8, alp_thresh, 1, 1, rho_BH,
                        eps_BH, vwlim_BH, sigma_max, inv_beta_max, true, false,
                        false, false, false, 1.0);
  c2p_1DPalenzuela c2p_Pal(eos_3p_ig, atmo, 100, 1e-8, alp_thresh, 1, 1, rho_BH,
                           eps_BH, vwlim_BH, sigma_max, inv_beta_max, true,
                           false, false, false, false, 1.0);
  c2p_1DRePrimAnd c2p_RPA(eos_3p_ig, atmo, 100, 1e-8, alp_thresh, 1, 1, rho_BH,
                          eps_BH, vwlim_BH, sigma_max, inv_beta_max, true,
                          false, false, false, false, 1.0);
  c2p_1DEntropy c2p_Ent(eos_3p_ig, atmo, 100, 1e-8, alp_thresh, 1, 1, rho_BH,
                        eps_BH, vwlim_BH, sigma_max, inv_beta_max, true, false,
                        false, false, false, 1.0);

  // Construct error report object:
  c2p_report rep_Noble;
  c2p_report rep_Pal;
  c2p_report rep_RPA;
  c2p_report rep_Ent;

  // Set primitive seeds
  const CCTK_REAL rho_in = 0.125;
  CCTK_REAL eps_in = 0.8;
  const CCTK_REAL Ye_in = 0.5;
  const CCTK_REAL press_in =
      eos_3p_ig->press_from_rho_eps_ye(rho_in, eps_in, Ye_in);
  const CCTK_REAL temp_in =
      eos_3p_ig->temp_from_rho_eps_ye(rho_in, eps_in, Ye_in);
  const CCTK_REAL entropy_in =
      eos_3p_ig->kappa_from_rho_eps_ye(rho_in, eps_in, Ye_in);
  const vec<CCTK_REAL, 3> vup_in = {0.0, 0.0, 0.0};
  const vec<CCTK_REAL, 3> Bup_in = {0.5, -0.5, 0.0};
  const vec<CCTK_REAL, 3> vdown_in = calc_contraction(g, vup_in);
  const CCTK_REAL wlor_in = calc_wlorentz(vdown_in, vup_in);

  prim_vars pv;
  // rho(p.I), eps(p.I), dummy_Ye, press(p.I), entropy, v_up, wlor, Bup
  prim_vars pv_seeds{rho_in,     eps_in, Ye_in,   press_in, temp_in,
                     entropy_in, vup_in, wlor_in, Bup_in};

  // cons_vars cv{dens(p.I), {momx(p.I), momy(p.I), momz(p.I)}, tau(p.I),
  //  dummy_DYe, DEnt, {dBx(p.I), dBy(p.I), dBz(p.I)}};

  cons_vars cv_Noble;
  cons_vars cv_Pal;
  cons_vars cv_Ent;
  cons_vars cv_all;

  // C2P solve may modify cons_vars input
  // Use the following to test C2Ps independently
  cv_Noble.from_prim(pv_seeds, g);
  cv_Pal.from_prim(pv_seeds, g);
  cv_Ent.from_prim(pv_seeds, g);

  // Use the following to test C2Ps in sequence
  // (As done in AsterX, the evolution thorn)
  cv_all.from_prim(pv_seeds, g);

  // Testing C2P Noble
  CCTK_VINFO("Testing C2P Noble...");
  c2p_Noble.solve(eos_3p_ig, pv, pv_seeds, cv_all, alp, beta, g, rep_Noble);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "temperature: %f, %f \n"
         "entropy: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.temperature, pv.temperature,
         pv_seeds.entropy, pv.entropy, pv_seeds.vel(0), pv.vel(0),
         pv_seeds.vel(1), pv.vel(1), pv_seeds.vel(2), pv.vel(2),
         pv_seeds.Bvec(0), pv.Bvec(0), pv_seeds.Bvec(1), pv.Bvec(1),
         pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "DYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n"
         "DEnt: %f \n",
         cv_all.dens, cv_all.tau, cv_all.mom(0), cv_all.mom(1), cv_all.mom(2),
         cv_all.DYe, cv_all.dBvec(0), cv_all.dBvec(1), cv_all.dBvec(2),
         cv_all.DEnt);
  /*
    assert(pv.rho == pv_seeds.rho);
    assert(pv.eps == pv_seeds.eps);
    assert(pv.press == pv_seeds.press);
    assert(pv.vel == pv_seeds.vel);
    assert(pv.Bvec == pv_seeds.Bvec);
  */

  rep_Noble.debug_message();

  // Testing C2P Palenzuela
  CCTK_VINFO("Testing C2P Palenzuela...");
  // c2p_Pal.solve(eos_3p_ig, pv, pv_seeds, cv, g, rep_Pal);
  c2p_Pal.solve(eos_3p_ig, pv, cv_all, alp, beta, g, rep_Pal);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "temperature: %f, %f \n"
         "entropy: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.temperature, pv.temperature,
         pv_seeds.entropy, pv.entropy, pv_seeds.vel(0), pv.vel(0),
         pv_seeds.vel(1), pv.vel(1), pv_seeds.vel(2), pv.vel(2),
         pv_seeds.Bvec(0), pv.Bvec(0), pv_seeds.Bvec(1), pv.Bvec(1),
         pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "DYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n"
         "DEnt: %f \n",
         cv_all.dens, cv_all.tau, cv_all.mom(0), cv_all.mom(1), cv_all.mom(2),
         cv_all.DYe, cv_all.dBvec(0), cv_all.dBvec(1), cv_all.dBvec(2),
         cv_all.DEnt);
  /*
    assert(pv.rho == pv_seeds.rho);
    assert(pv.eps == pv_seeds.eps);
    assert(pv.press == pv_seeds.press);
    assert(pv.vel == pv_seeds.vel);
    assert(pv.Bvec == pv_seeds.Bvec);
  */
  rep_Pal.debug_message();

  // Testing C2P RePrimAnd
  CCTK_VINFO("Testing C2P RePrimAnd...");
  c2p_RPA.solve(eos_3p_ig, pv, cv_all, alp, beta, g, rep_RPA);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "temperature: %f, %f \n"
         "entropy: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.temperature, pv.temperature,
         pv_seeds.entropy, pv.entropy, pv_seeds.vel(0), pv.vel(0),
         pv_seeds.vel(1), pv.vel(1), pv_seeds.vel(2), pv.vel(2),
         pv_seeds.Bvec(0), pv.Bvec(0), pv_seeds.Bvec(1), pv.Bvec(1),
         pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "DYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n"
         "DEnt: %f \n",
         cv_all.dens, cv_all.tau, cv_all.mom(0), cv_all.mom(1), cv_all.mom(2),
         cv_all.DYe, cv_all.dBvec(0), cv_all.dBvec(1), cv_all.dBvec(2),
         cv_all.DEnt);
  rep_RPA.debug_message();

  // Testing C2P Entropy
  CCTK_VINFO("Testing C2P Entropy...");
  c2p_Ent.solve(eos_3p_ig, pv, cv_all, alp, beta, g, rep_Ent);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "temperature: %f, %f \n"
         "entropy: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.temperature, pv.temperature,
         pv_seeds.entropy, pv.entropy, pv_seeds.vel(0), pv.vel(0),
         pv_seeds.vel(1), pv.vel(1), pv_seeds.vel(2), pv.vel(2),
         pv_seeds.Bvec(0), pv.Bvec(0), pv_seeds.Bvec(1), pv.Bvec(1),
         pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "DYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n"
         "DEnt: %f \n",
         cv_all.dens, cv_all.tau, cv_all.mom(0), cv_all.mom(1), cv_all.mom(2),
         cv_all.DYe, cv_all.dBvec(0), cv_all.dBvec(1), cv_all.dBvec(2),
         cv_all.DEnt);
  /*
    assert(pv.rho == pv_seeds.rho);
    assert(pv.eps == pv_seeds.eps);
    assert(pv.press == pv_seeds.press);
    assert(pv.vel == pv_seeds.vel);
    assert(pv.Bvec == pv_seeds.Bvec);
  */
  rep_Ent.debug_message();
}

} // namespace Con2PrimFactory
