#ifndef C2P_1DPALENZUELA_HXX
#define C2P_1DPALENZUELA_HXX

#include "c2p.hxx"
#include "roots.hxx"
namespace Con2PrimFactory {

class c2p_1DPalenzuela : public c2p {
public:
  /* Some attributes */
  CCTK_REAL GammaIdealFluid;

  /* Constructor */
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DPalenzuela(
      EOSType &eos_th, atmosphere &atm, CCTK_INT maxIter, CCTK_REAL tol,
      CCTK_REAL rho_str, CCTK_REAL vwlim, CCTK_REAL B_lim, bool ye_len);

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
      get_Ssq_Exact(vec<CCTK_REAL, 3> &mom,
                    const smat<CCTK_REAL, 3> &gup) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
      get_Bsq_Exact(vec<CCTK_REAL, 3> &B_up,
                    const smat<CCTK_REAL, 3> &glo) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
      get_BiSi_Exact(vec<CCTK_REAL, 3> &Bvec, vec<CCTK_REAL, 3> &mom) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 3>
      get_WLorentz_vsq_bsq_Seeds(vec<CCTK_REAL, 3> &B_up,
                                 vec<CCTK_REAL, 3> &v_up,
                                 const smat<CCTK_REAL, 3> &glo) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_to_nan(prim_vars &pv, cons_vars &cv) const;
  /* Called by 1DPalenzuela */
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  xPalenzuelaToPrim(CCTK_REAL xPalenzuela_Sol, CCTK_REAL Ssq, CCTK_REAL Bsq,
                    CCTK_REAL BiSi, EOSType &eos_th, prim_vars &pv,
                    cons_vars cv, const smat<CCTK_REAL, 3> &gup,
                    const smat<CCTK_REAL, 3> &glo) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  funcRoot_1DPalenzuela(CCTK_REAL Ssq, CCTK_REAL Bsq, CCTK_REAL BiSi,
                        CCTK_REAL x, EOSType &eos_th, cons_vars &cv) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  solve(EOSType &eos_th, prim_vars &pv, prim_vars &pv_seeds, cons_vars cv,
        const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const;

  /* Destructor */
  CCTK_HOST CCTK_DEVICE ~c2p_1DPalenzuela();
};

/* Constructor */
template <typename EOSType>
CCTK_HOST CCTK_DEVICE
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DPalenzuela::c2p_1DPalenzuela(
        EOSType &eos_th, atmosphere &atm, CCTK_INT maxIter, CCTK_REAL tol,
        CCTK_REAL rho_str, CCTK_REAL vwlim, CCTK_REAL B_lim, bool ye_len) {

  GammaIdealFluid = eos_th.gamma;
  maxIterations = maxIter;
  tolerance = tol;
  rho_strict = rho_str;
  ye_lenient = ye_len;
  vw_lim = vwlim;
  w_lim = sqrt(1.0 + vw_lim * vw_lim);
  v_lim = vw_lim / w_lim;
  Bsq_lim = B_lim * B_lim;
  atmo = atm;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DPalenzuela::get_Ssq_Exact(vec<CCTK_REAL, 3> &mom,
                                const smat<CCTK_REAL, 3> &gup) const {

  CCTK_REAL Ssq;
  /* calculate S_squared */
  Ssq = mom(X) * (gup(X, X) * mom(X) + gup(X, Y) * mom(Y) + gup(X, Z) * mom(Z));
  Ssq +=
      mom(Y) * (gup(X, Y) * mom(X) + gup(Y, Y) * mom(Y) + gup(Y, Z) * mom(Z));
  Ssq +=
      mom(Z) * (gup(X, Z) * mom(X) + gup(Y, Z) * mom(Y) + gup(Z, Z) * mom(Z));

  if ((Ssq < 0.) && (fabs(Ssq) < 1.0e-13)) {
    Ssq = fabs(Ssq);
  }
  return Ssq;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DPalenzuela::get_Bsq_Exact(vec<CCTK_REAL, 3> &B_up,
                                const smat<CCTK_REAL, 3> &glo) const {

  vec<CCTK_REAL, 3> B_low = calc_contraction(glo, B_up);
  return calc_contraction(B_low, B_up); // Bsq
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DPalenzuela::get_BiSi_Exact(vec<CCTK_REAL, 3> &Bvec,
                                 vec<CCTK_REAL, 3> &mom) const {

  return Bvec(X) * mom(X) + Bvec(Y) * mom(Y) + Bvec(Z) * mom(Z); // BiSi
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 3>
c2p_1DPalenzuela::get_WLorentz_vsq_bsq_Seeds(
    vec<CCTK_REAL, 3> &B_up, vec<CCTK_REAL, 3> &v_up,
    const smat<CCTK_REAL, 3> &glo) const {
  vec<CCTK_REAL, 3> v_low = calc_contraction(glo, v_up);
  CCTK_REAL vsq = calc_contraction(v_low, v_up);
  CCTK_REAL VdotB = calc_contraction(v_low, B_up);
  CCTK_REAL VdotBsq = VdotB * VdotB;
  CCTK_REAL Bsq = get_Bsq_Exact(B_up, glo);

  if ((vsq < 0.) && (fabs(vsq) < 1.0e-13)) {
    vsq = fabs(vsq);
  }

  CCTK_REAL w_lor = 1. / sqrt(1. - vsq);
  CCTK_REAL bsq = ((Bsq) / (w_lor * w_lor)) + VdotBsq;
  vec<CCTK_REAL, 3> w_vsq_bsq{w_lor, vsq, bsq};

  return w_vsq_bsq; //{w_lor, vsq, bsq}
}

/* Called by 1DPalenzuela */

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_1DPalenzuela::xPalenzuelaToPrim(CCTK_REAL xPalenzuela_Sol, CCTK_REAL Ssq,
                                    CCTK_REAL Bsq, CCTK_REAL BiSi,
                                    EOSType &eos_th, prim_vars &pv,
                                    cons_vars cv, const smat<CCTK_REAL, 3> &gup,
                                    const smat<CCTK_REAL, 3> &glo) const {
  const CCTK_REAL qPalenzuela = cv.tau / cv.dens;
  const CCTK_REAL rPalenzuela = Ssq / pow(cv.dens, 2);
  const CCTK_REAL sPalenzuela = Bsq / cv.dens;
  const CCTK_REAL tPalenzuela = BiSi / pow(cv.dens, 3. / 2.);

  // (i)
  CCTK_REAL Wminus2 =
      1.0 -
      (xPalenzuela_Sol * xPalenzuela_Sol * rPalenzuela +
       (2.0 * xPalenzuela_Sol + sPalenzuela) * tPalenzuela * tPalenzuela) /
          (xPalenzuela_Sol * xPalenzuela_Sol * (xPalenzuela_Sol + sPalenzuela) *
           (xPalenzuela_Sol + sPalenzuela));
  Wminus2 = fmin(fmax(Wminus2, 1e-10), 1 - 1e-10);
  const CCTK_REAL W_sol = pow(Wminus2, -0.5);

  // (ii)
  pv.rho = cv.dens / W_sol;

  // (iii)
  pv.eps = W_sol - 1.0 + (1.0 - W_sol * W_sol) * xPalenzuela_Sol / W_sol +
           W_sol * (qPalenzuela - sPalenzuela +
                    tPalenzuela * tPalenzuela /
                        (2 * xPalenzuela_Sol * xPalenzuela_Sol) +
                    sPalenzuela / (2.0 * W_sol * W_sol));

  // (iv)
  // CCTK_REAL P_loc = get_Press_funcRhoEps(rho_loc, eps_loc);

  // Taken from WZ2Prim (2DNRNoble)
  CCTK_REAL Z_Sol = xPalenzuela_Sol * pv.rho * W_sol;

  pv.vel(X) =
      (gup(X, X) * cv.mom(X) + gup(X, Y) * cv.mom(Y) + gup(X, Z) * cv.mom(Z)) /
      (Z_Sol + Bsq);
  pv.vel(X) += BiSi * cv.dBvec(X) / (Z_Sol * (Z_Sol + Bsq));

  pv.vel(Y) =
      (gup(X, Y) * cv.mom(X) + gup(Y, Y) * cv.mom(Y) + gup(Y, Z) * cv.mom(Z)) /
      (Z_Sol + Bsq);
  pv.vel(Y) += BiSi * cv.dBvec(Y) / (Z_Sol * (Z_Sol + Bsq));

  pv.vel(Z) =
      (gup(X, Z) * cv.mom(X) + gup(Y, Z) * cv.mom(Y) + gup(Z, Z) * cv.mom(Z)) /
      (Z_Sol + Bsq);
  pv.vel(Z) += BiSi * cv.dBvec(Z) / (Z_Sol * (Z_Sol + Bsq));

  pv.w_lor = W_sol;

  pv.Ye = cv.dYe / cv.dens;

  pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

  pv.Bvec = cv.dBvec;

  const vec<CCTK_REAL, 3> Elow = calc_cross_product(pv.Bvec, pv.vel);
  pv.E = calc_contraction(gup, Elow);
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DPalenzuela::funcRoot_1DPalenzuela(CCTK_REAL Ssq, CCTK_REAL Bsq,
                                        CCTK_REAL BiSi, CCTK_REAL x,
                                        EOSType &eos_th, cons_vars &cv) const {
  // computes f(x) from x and q,r,s,t
  const CCTK_REAL qPalenzuela = cv.tau / cv.dens;
  const CCTK_REAL rPalenzuela = Ssq / pow(cv.dens, 2);
  const CCTK_REAL sPalenzuela = Bsq / cv.dens;
  const CCTK_REAL tPalenzuela = BiSi / pow(cv.dens, 3. / 2.);

  // (i)
  CCTK_REAL Wminus2 =
      1.0 - (x * x * rPalenzuela +
             (2.0 * x + sPalenzuela) * tPalenzuela * tPalenzuela) /
                (x * x * (x + sPalenzuela) * (x + sPalenzuela));
  Wminus2 = fmin(fmax(Wminus2, 1e-10), 1 - 1e-10);
  const CCTK_REAL W_loc = pow(Wminus2, -0.5);

  // (ii)
  CCTK_REAL rho_loc = cv.dens / W_loc;
  CCTK_REAL Ye_loc = cv.dYe / cv.dens;

  // (iii)
  CCTK_REAL eps_loc = W_loc - 1.0 + (1.0 - W_loc * W_loc) * x / W_loc +
                      W_loc * (qPalenzuela - sPalenzuela +
                               tPalenzuela * tPalenzuela / (2 * x * x) +
                               sPalenzuela / (2 * W_loc * W_loc));

  // (iv)
  CCTK_REAL P_loc =
      eos_th.press_from_valid_rho_eps_ye(rho_loc, eps_loc, Ye_loc);

  return (x - (1.0 + eps_loc + P_loc / rho_loc) * W_loc);
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_1DPalenzuela::solve(EOSType &eos_th, prim_vars &pv, prim_vars &pv_seeds,
                        cons_vars cv, const smat<CCTK_REAL, 3> &glo,
                        c2p_report &rep) const {

  ROOTSTAT status = ROOTSTAT::SUCCESS;
  rep.iters = 0;
  rep.adjust_cons = false;
  rep.set_atmo = false;
  rep.status = c2p_report::SUCCESS;

  /* Check validity of the 3-metric and compute its inverse */
  const CCTK_REAL spatial_detg = calc_det(glo);
  const CCTK_REAL sqrt_detg = sqrt(spatial_detg);
  if ((!isfinite(sqrt_detg)) || (sqrt_detg <= 0)) {
    rep.set_invalid_detg(sqrt_detg);
    set_to_nan(pv, cv);
    return;
  }
  const smat<CCTK_REAL, 3> gup = calc_inv(glo, spatial_detg);

  /* Undensitize the conserved vars */
  cv.dens /= sqrt_detg;
  cv.tau /= sqrt_detg;
  cv.mom /= sqrt_detg;
  cv.dBvec /= sqrt_detg;
  cv.dYe /= sqrt_detg;

  if (cv.dens <= atmo.rho_cut) {
    rep.set_atmo_set();
    atmo.set(pv, cv, glo);
    return;
  }

  // compute primitive B seed from conserved B of current time step for better
  // guess
  pv_seeds.Bvec = cv.dBvec;

  const CCTK_REAL Ssq = get_Ssq_Exact(cv.mom, gup);
  const CCTK_REAL Bsq = get_Bsq_Exact(pv_seeds.Bvec, glo);
  const CCTK_REAL BiSi = get_BiSi_Exact(pv_seeds.Bvec, cv.mom);
  const vec<CCTK_REAL, 3> w_vsq_bsq = get_WLorentz_vsq_bsq_Seeds(
      pv_seeds.Bvec, pv_seeds.vel, glo); // this also recomputes pv_seeds.w_lor
  pv_seeds.w_lor = w_vsq_bsq(0);
  // CCTK_REAL vsq_seed = w_vsq_bsq(1);
  // const CCTK_REAL bsq = w_vsq_bsq(2);

  if ((!isfinite(cv.dens)) || (!isfinite(Ssq)) || (!isfinite(Bsq)) ||
      (!isfinite(BiSi)) || (!isfinite(cv.dYe))) {
    rep.set_nans_in_cons(cv.dens, Ssq, Bsq, BiSi, cv.dYe);
    set_to_nan(pv, cv);
    return;
  }

  if (Bsq < 0) {
    rep.set_neg_Bsq(Bsq);
    set_to_nan(pv, cv);
    return;
  }

  if (Bsq > Bsq_lim) {
    rep.set_B_limit(Bsq);
    set_to_nan(pv, cv);
    return;
  }

  /* update rho seed from cv and wlor */
  // rho consistent with cv.rho should be better guess than rho from last
  // timestep
  pv_seeds.rho = cv.dens / pv_seeds.w_lor;

  // Find x, this is the recovery process
  const CCTK_INT minbits = std::numeric_limits<CCTK_REAL>::digits - 4;
  const CCTK_INT maxiters = maxIterations;

  CCTK_REAL qPalenzuela = cv.tau / cv.dens;
  CCTK_REAL sPalenzuela = Bsq / cv.dens;
  CCTK_REAL xPalenzuela_lowerBound = 1.0 + qPalenzuela - sPalenzuela;
  CCTK_REAL xPalenzuela_upperBound = 2.0 + 2.0 * qPalenzuela - sPalenzuela;
  CCTK_REAL a = xPalenzuela_lowerBound;
  CCTK_REAL b = xPalenzuela_upperBound;
  auto fn = [&](auto x) {
    return funcRoot_1DPalenzuela(Ssq, Bsq, BiSi, x, eos_th, cv);
  };
  auto result = Algo::brent(fn, a, b, minbits, maxiters, rep.iters);

  // Pick best solution
  CCTK_REAL xPalenzuela_Sol;
  if (abs(fn(result.first)) < abs(fn(result.second))) {
    xPalenzuela_Sol = result.first;
  } else {
    xPalenzuela_Sol = result.second;
  }

  // Check solution and calculate primitives
  // TODO:check if to pass result.first or xPalenzuela_Sol
  if (rep.iters < maxiters && abs(fn(xPalenzuela_Sol)) < tolerance) {
    rep.status = c2p_report::SUCCESS;
    status = ROOTSTAT::SUCCESS;
  } else {
    // set status to root not converged
    rep.set_root_conv();
    status = ROOTSTAT::NOT_CONVERGED;
  }

  xPalenzuelaToPrim(xPalenzuela_Sol, Ssq, Bsq, BiSi, eos_th, pv, cv, gup, glo);

  // set to atmo if computed rho is below floor density
  if (pv.rho < atmo.rho_cut) {
    rep.set_atmo_set();
    atmo.set(pv, cv, glo);
    return;
  }

  // check the validity of the computed eps
  auto rgeps = eos_th.range_eps_from_valid_rho_ye(pv.rho, pv.Ye);
  if (pv.eps > rgeps.max) {
    //printf("(pv.eps > rgeps.max) is true, adjusting cons.. \n");
    rep.adjust_cons = true;
    if (pv.rho >= rho_strict) {
      rep.set_range_eps(pv.eps); // sets adjust_cons to false by default
      rep.adjust_cons = true;
      set_to_nan(pv, cv);
      return;
    }
  } else if (pv.eps < rgeps.min) {
    /*
    printf(
        "(pv.eps < rgeps.min) is true! pv.eps, rgeps.min: %26.16e, %26.16e \n",
        pv.eps, rgeps.min);
    printf(" Not adjusting cons.. \n");
    */
    rep.set_range_eps(rgeps.min); // sets adjust_cons to true
  }

  // TODO: check validity for Ye

  // check if computed velocities are within the specified limit
  vec<CCTK_REAL, 3> v_low = calc_contraction(glo, pv.vel);
  CCTK_REAL vsq_Sol = calc_contraction(v_low, pv.vel);
  CCTK_REAL sol_v = sqrt(vsq_Sol);
  if (sol_v > v_lim) {
    /*
    printf("(sol_v > v_lim) is true! \n");
    printf("sol_v, v_lim: %26.16e, %26.16e \n", sol_v, v_lim);
    */
    pv.rho = cv.dens / w_lim;
    if (pv.rho >= rho_strict) {
      rep.set_speed_limit({sol_v, sol_v, sol_v});
      set_to_nan(pv, cv);
      return;
    }
    pv.vel *= v_lim / sol_v;
    pv.w_lor = w_lim;
    pv.eps = std::min(std::max(eos_th.rgeps.min, pv.eps), eos_th.rgeps.max);
    pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  }

  // Recompute cons if prims have been adjusted
  if (rep.adjust_cons) {
    cv.from_prim(pv, glo);
  }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_1DPalenzuela::set_to_nan(prim_vars &pv, cons_vars &cv) const {
  pv.set_to_nan();
  cv.set_to_nan();
}

/* Destructor */
CCTK_HOST CCTK_DEVICE
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DPalenzuela::~c2p_1DPalenzuela() {
  // How to destruct properly a vector?
}
} // namespace Con2PrimFactory

#endif
