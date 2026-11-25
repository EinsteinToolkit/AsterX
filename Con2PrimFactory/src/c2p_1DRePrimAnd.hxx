#ifndef C2P_1DREPRIMAND_HXX
#define C2P_1DREPRIMAND_HXX

#include <algorithm>
#include <cmath>

#include "c2p.hxx"
#include "c2p_report.hxx"
#include "prims.hxx"
#include "cons.hxx"
#include "roots.hxx"

#include "c2p_1DRePrimAnd_intervals.hxx"
#include "c2p_1DRePrimAnd_rootfinder.hxx"
#include "c2p_1DRePrimAnd_internals.hxx"

namespace Con2PrimFactory {

class c2p_1DRePrimAnd : public c2p {
public:

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DRePrimAnd(
      const EOSType *eos_3p, const atmosphere &atm_in, CCTK_INT maxIter,
      CCTK_REAL tol, CCTK_REAL alp_thresh_in, CCTK_REAL vwlim,
      CCTK_REAL B_lim, CCTK_REAL rho_BH_in, CCTK_REAL eps_BH_in,
      CCTK_REAL vwlim_BH_in, CCTK_REAL sigma_max_in,
      CCTK_REAL inv_beta_max_in, bool ye_len, bool use_z,
      bool use_temperature, bool use_pressure_atmo, bool soft_root_conv,
      CCTK_REAL soft_root_width_factor_in) {

    atmo = atm_in;
    maxIterations = maxIter;
    tolerance = tol;
    alp_thresh = alp_thresh_in;
    vw_lim = vwlim;
    w_lim = sqrt(1.0 + vw_lim * vw_lim);
    v_lim = vw_lim / w_lim;
    Bsq_lim = B_lim * B_lim;
    rho_BH = rho_BH_in;
    eps_BH = eps_BH_in;
    vwlim_BH = vwlim_BH_in;
    ye_lenient = ye_len;
    use_zprim = use_z;
    use_temp = use_temperature;
    use_press_atmo = use_pressure_atmo;
    soft_root_convergence = soft_root_conv;
    soft_root_width_factor = fmax(CCTK_REAL(1.0), soft_root_width_factor_in);
    sigma_max = sigma_max_in;
    inv_beta_max = inv_beta_max_in;

  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_to_nan(prim_vars &pv, cons_vars &cv) const {
    pv.set_to_nan();
    cv.set_to_nan();
  }

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  solve(const EOSType *eos_3p, prim_vars &pv, cons_vars &cv,
        const CCTK_REAL alp, const vec<CCTK_REAL, 3> &beta,
        const smat<CCTK_REAL, 3> &glo, c2p_report &rep,
        bool reject_nonpositive_eps = false) const {
    rep.iters = 0;
    rep.adjust_cons = false;
    rep.set_atmo = false;
    rep.soft_root_conv = false;
    rep.status = c2p_report::SUCCESS;

    const CCTK_REAL spatial_detg = calc_det(glo);
    const CCTK_REAL sqrt_detg = sqrt(spatial_detg);
    const bool minor1{glo(X, X) > 0.0};
    const bool minor2{glo(X, X) * glo(Y, Y) - glo(X, Y) * glo(X, Y) > 0.0};
    const bool minor3{spatial_detg > 0.0};
    if (!(minor1 && minor2 && minor3)) {
      rep.set_invalid_detg(sqrt_detg);
      set_to_nan(pv, cv);
      return;
    }

    const smat<CCTK_REAL, 3> gup = calc_inv(glo, spatial_detg);

    const cons_vars cv_const = cv;

    cv.dens /= sqrt_detg;
    cv.tau /= sqrt_detg;
    cv.mom /= sqrt_detg;
    cv.dBvec /= sqrt_detg;
    cv.DYe /= sqrt_detg;
    cv.DEnt /= sqrt_detg;

    const CCTK_REAL Ssq =
        calc_contraction(calc_contraction(gup, cv.mom), cv.mom);
    const CCTK_REAL Bsq =
        calc_contraction(calc_contraction(glo, cv.dBvec), cv.dBvec);
    const CCTK_REAL BiSi = cv.dBvec(X) * cv.mom(X) + cv.dBvec(Y) * cv.mom(Y) +
                           cv.dBvec(Z) * cv.mom(Z);

    if ((!isfinite(cv.dens)) || (!isfinite(Ssq)) || (!isfinite(Bsq)) ||
        (!isfinite(BiSi)) || (!isfinite(cv.DYe)) || (!isfinite(cv.DEnt))) {
      rep.set_nans_in_cons(cv.dens, Ssq, Bsq, BiSi, cv.DYe);
      set_to_nan(pv, cv);
      return;
    }

    if (Bsq > Bsq_lim) {
      rep.set_B_limit(Bsq);
      set_to_nan(pv, cv);
      return;
    }

    // Ye needs to be a non-const lvalue for EOS calls
    CCTK_REAL Ye_raw = cv.DYe / cv.dens;
    const CCTK_REAL Ye = fmin(fmax(eos_3p->rgye.min, Ye_raw), eos_3p->rgye.max);
    const bool ye_clipped =
        (Ye_raw < eos_3p->rgye.min) || (Ye_raw > eos_3p->rgye.max);

    // A positive enthalpy lower bound requires eps_min > -1.
    if ((!isfinite(eos_3p->rgeps.min)) || eos_3p->rgeps.min <= -1.0) {
      rep.set_range_eps(eos_3p->rgeps.min);
      cv = cv_const;
      return;
    }

    typename RePrimAnd::froot<EOSType>::cache cache{};
    RePrimAnd::froot<EOSType> f(
        eos_3p, Ye, cv.dens, cv.tau / cv.dens, Ssq / (cv.dens * cv.dens),
        (BiSi * BiSi) / (cv.dens * cv.dens * cv.dens), Bsq / cv.dens, cache);

    ROOTSTAT prep_stat{};
    interval<CCTK_REAL> mu_br = f.initial_bracket(prep_stat);
    if (prep_stat != ROOTSTAT::SUCCESS) {
      if (prep_stat == ROOTSTAT::NOT_BRACKETED)
        rep.set_prep_root_bracket();
      else
        rep.set_prep_root_conv();
      cv = cv_const;
      return;
    }

    interval<CCTK_REAL> rgrho_iv{eos_3p->rgrho.min, eos_3p->rgrho.max};
    RePrimAnd::rarecase<EOSType> rc(mu_br, rgrho_iv, f);
    if (rc.rho_too_big) {
      rep.set_range_rho(cv.dens, cv.dens);
      cv = cv_const;
      return;
    }
    if (rc.rho_too_small) {
      rep.set_atmo_set();
      pv.Bvec = cv.dBvec;
      atmo.set(pv, cv, glo);
      return;
    }
    mu_br = rc.bracket;

    const CCTK_INT maxiters = maxIterations;

    ROOTSTAT root_stat{};
    interval<CCTK_REAL> result =
        findroot_no_deriv(f, mu_br, tolerance, (unsigned int)maxiters, root_stat);
    rep.iters = (CCTK_INT)cache.calls;

    if (root_stat == ROOTSTAT::NOT_CONVERGED) {
      bool accept_soft = false;
      if (soft_root_convergence) {
        const CCTK_REAL a = fmin(result.min(), result.max());
        const CCTK_REAL b = fmax(result.min(), result.max());
        if (std::isfinite(a) && std::isfinite(b) && b >= a) {
          const CCTK_REAL width = b - a;
          const CCTK_REAL scale =
              fmax(CCTK_REAL(0.0), fmax(std::abs(a), std::abs(b)));
          const CCTK_REAL width_tol =
              soft_root_width_factor * tolerance * scale;
          const CCTK_REAL mu_ref =
              fmax(CCTK_REAL(0.5) * (a + b), CCTK_REAL(1.0e-300));
          const bool stopif_ok = f.stopif(mu_ref, width, tolerance);
          accept_soft = (width <= width_tol) || stopif_ok;
        }
      }
      if (!accept_soft) {
        rep.set_root_conv();
        cv = cv_const;
        return;
      }
      rep.set_soft_root_conv();
      root_stat = ROOTSTAT::SUCCESS;
    }

    if (root_stat != ROOTSTAT::SUCCESS) {
      if (rc.rho_big) {
        rep.set_range_rho(cv.dens, cv.dens);
        cv = cv_const;
        return;
      }
      if (rc.rho_small) {
        rep.set_atmo_set();
        pv.Bvec = cv.dBvec;
        atmo.set(pv, cv, glo);
        return;
      }
      rep.set_root_bracket();
      cv = cv_const;
      return;
    }

    CCTK_REAL mu = cache.lmu;
    const CCTK_REAL mu_lo = fmin(result.min(), result.max());
    const CCTK_REAL mu_hi = fmax(result.min(), result.max());

    if (!std::isfinite(mu) || mu < mu_lo || mu > mu_hi) {
      auto fn = [&](CCTK_REAL x) { return f(x); };
      const CCTK_REAL a_root = result.min();
      const CCTK_REAL b_root = result.max();
      const CCTK_REAL fa = fn(a_root);
      const CCTK_REAL fb = fn(b_root);
      mu = (fb == (CCTK_REAL)0 || std::abs(fb) < std::abs(fa)) ? b_root
           : (std::abs(fa) < std::abs(fb))                     ? a_root
                                     : CCTK_REAL(0.5) * (a_root + b_root);

      // Ensure cache corresponds to the final chosen mu.
      (void)f(mu);
    }

    if ((!isfinite(cache.rho)) || cache.rho <= 0.0) {
      rep.set_range_rho(cv.dens, cache.rho);
      cv = cv_const;
      return;
    }

    const auto rgeps = eos_3p->range_eps_from_rho_ye(cache.rho, cache.ye);
    const bool eps_clipped =
        cache.eps_raw < rgeps.min || cache.eps_raw > rgeps.max;

    // Negative energy is valid only within the local EOS energy range.
    const bool eps_invalid =
        rgeps.min < 0.0 ? cache.eps_raw < rgeps.min : cache.eps_raw <= 0.0;
    if ((!isfinite(cache.eps_raw)) ||
        (reject_nonpositive_eps && eps_invalid)) {
      rep.set_range_eps(cache.eps_raw);
      cv = cv_const;
      return;
    }

    // ------------------------------------------------------------------
    // Use the EOS-consistent RePrimAnd cached primitives
    // ------------------------------------------------------------------
    pv.rho = cache.rho;
    pv.Ye = cache.ye;
    pv.eps = cache.eps;
    pv.press = cache.press;
    pv.w_lor = cache.w;

    // If root-cache thermodynamics were clipped to EOS bounds, force
    // conservative recomputation for consistency with adjusted primitives.
    if (ye_clipped || eps_clipped) {
      rep.adjust_cons = true;
    }

    pv.temperature = eos_3p->temp_from_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.entropy = eos_3p->kappa_from_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    // ------------------------------------------------------------------
    // Velocity: RePrimAnd magnetic-aware reconstruction
    //   v^i = mu * x * ( r^i + (r.b) * mu * b^i )
    // where r^i = S^i / D, b^i = B^i / sqrt(D), (r.b) = (B^i S_i)/(D^(3/2)).
    // ------------------------------------------------------------------

    const vec<CCTK_REAL, 3> mom_up = calc_contraction(gup, cv.mom);

    const CCTK_REAL D = cv.dens;
    const CCTK_REAL sqD = std::sqrt(D);

    const vec<CCTK_REAL, 3> r_u = mom_up * (1.0 / (D + 1e-300));
    const vec<CCTK_REAL, 3> b_u = cv.dBvec * (1.0 / (sqD + 1e-300));
    const CCTK_REAL rb = BiSi / (D * (sqD + 1e-300));

    pv.vel = (mu * cache.x) * (r_u + (rb * mu) * b_u);

    const CCTK_REAL sol_v = std::sqrt(std::max(CCTK_REAL(0), cache.vsqr));
    if (sol_v > v_lim) {
      pv.rho = cv.dens / w_lim;
      pv.vel *= v_lim / (sol_v + CCTK_REAL(1e-300));
      pv.w_lor = w_lim;
      const auto rgeps_lim = eos_3p->range_eps_from_rho_ye(pv.rho, pv.Ye);
      pv.eps = fmin(fmax(rgeps_lim.min, pv.eps), rgeps_lim.max);
      pv.press = eos_3p->press_from_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
      pv.temperature = eos_3p->temp_from_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
      pv.entropy = eos_3p->kappa_from_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
      rep.adjust_cons = true;
    } else {
      pv.w_lor = cache.w;
    }

    pv.Bvec = cv.dBvec;
    const vec<CCTK_REAL, 3> Elow = calc_cross_product(pv.Bvec, pv.vel);
    pv.E = calc_contraction(gup, Elow);

    // set to atmo if computed rho is below floor density
    // and atmo obeys magnetic field limits
    const CCTK_REAL b2_atm = calc_norm(pv.Bvec, glo);
    if (pv.rho < atmo.rho_cut &&
        (b2_atm / atmo.rho_atmo <= sigma_max &&
         b2_atm / (2.0 * atmo.press_atmo) <= inv_beta_max)) {
      rep.set_atmo_set();
      atmo.set(pv, cv, glo);
      return;
    }

    c2p::prims_floors_and_ceilings(eos_3p, pv, cv, alp, beta, glo, rep);
    if (rep.failed()) {
      cv = cv_const;
      return;
    }

    if (rep.adjust_cons) {
      cv.from_prim(pv, glo);
      cv.dBvec = cv_const.dBvec;
    } else {
      cv = cv_const;
      cv.DEnt = cv.dens * pv.entropy;
    }
  }

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  solve(const EOSType *eos_3p, prim_vars &pv, cons_vars &cv,
        const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const {
    // Backward-compatible overload for legacy call sites.
    const CCTK_REAL alp = 1.0;
    const vec<CCTK_REAL, 3> beta{0.0, 0.0, 0.0};
    solve(eos_3p, pv, cv, alp, beta, glo, rep);
  }

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  operator()(const EOSType *eos_3p, prim_vars &pv, cons_vars &cv,
             const smat<CCTK_REAL, 3> & /*gup_unused*/,
             const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const {
    solve(eos_3p, pv, cv, glo, rep);
  }
};

} // namespace Con2PrimFactory

#endif
