#ifndef C2P_1DENTROPY_HXX
#define C2P_1DENTROPY_HXX

#include "c2p.hxx"
#include "roots.hxx"
namespace Con2PrimFactory {

class c2p_1DEntropy : public c2p {
public:
  /* Some attributes */
  CCTK_REAL GammaIdealFluid;

  /* Constructor */
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DEntropy(
      const EOSType *eos_3p, const atmosphere &atm, CCTK_INT maxIter,
      CCTK_REAL tol, CCTK_REAL alp_thresh_in, CCTK_REAL consError,
      CCTK_REAL vwlim, CCTK_REAL B_lim, CCTK_REAL rho_BH_in,
      CCTK_REAL eps_BH_in, CCTK_REAL vwlim_BH_in, bool ye_len, bool use_z,
      bool use_temperature);

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Ssq_Exact(const vec<CCTK_REAL, 3> &mom,
                const smat<CCTK_REAL, 3> &gup) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Bsq_Exact(const vec<CCTK_REAL, 3> &B_up,
                const smat<CCTK_REAL, 3> &glo) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_BiSi_Exact(const vec<CCTK_REAL, 3> &Bvec,
                 const vec<CCTK_REAL, 3> &mom) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 3>
  get_WLorentz_vsq_bsq_Seeds(const vec<CCTK_REAL, 3> &B_up,
                             const vec<CCTK_REAL, 3> &v_up,
                             const smat<CCTK_REAL, 3> &glo) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_to_nan(prim_vars &pv, cons_vars &cv) const;
  /* Called by 1DEntropy */
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  xEntropyToPrim(CCTK_REAL xEntropy_Sol, CCTK_REAL Ssq, CCTK_REAL Bsq,
                 CCTK_REAL BiSi, const EOSType *eos_3p, prim_vars &pv,
                 const cons_vars &cv, const smat<CCTK_REAL, 3> &gup,
                 const smat<CCTK_REAL, 3> &glo) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  funcRoot_1DEntropy(CCTK_REAL Ssq, CCTK_REAL Bsq, CCTK_REAL BiSi, CCTK_REAL x,
                     const EOSType *eos_3p, const cons_vars &cv) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  solve(const EOSType *eos_3p, prim_vars &pv, cons_vars &cv,
        const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const;

  /* Destructor */
  CCTK_HOST CCTK_DEVICE ~c2p_1DEntropy();
};

/* Constructor */
template <typename EOSType>
CCTK_HOST CCTK_DEVICE
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DEntropy::c2p_1DEntropy(
        const EOSType *eos_3p, const atmosphere &atm, CCTK_INT maxIter,
        CCTK_REAL tol, CCTK_REAL alp_thresh_in, CCTK_REAL consError,
        CCTK_REAL vwlim, CCTK_REAL B_lim, CCTK_REAL rho_BH_in,
        CCTK_REAL eps_BH_in, CCTK_REAL vwlim_BH_in, bool ye_len, bool use_z,
        bool use_temperature) {

  // Base
  atmo = atm;
  maxIterations = maxIter;
  tolerance = tol;
  alp_thresh = alp_thresh_in;
  cons_error = consError;
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

  // Derived
  GammaIdealFluid = eos_3p->gamma;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DEntropy::get_Ssq_Exact(const vec<CCTK_REAL, 3> &mom,
                             const smat<CCTK_REAL, 3> &gup) const {

  CCTK_REAL Ssq;
  /* calculate S_squared */
  Ssq = mom(X) * (gup(X, X) * mom(X) + gup(X, Y) * mom(Y) + gup(X, Z) * mom(Z));
  Ssq +=
      mom(Y) * (gup(X, Y) * mom(X) + gup(Y, Y) * mom(Y) + gup(Y, Z) * mom(Z));
  Ssq +=
      mom(Z) * (gup(X, Z) * mom(X) + gup(Y, Z) * mom(Y) + gup(Z, Z) * mom(Z));

  return Ssq;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DEntropy::get_Bsq_Exact(const vec<CCTK_REAL, 3> &B_up,
                             const smat<CCTK_REAL, 3> &glo) const {

  vec<CCTK_REAL, 3> B_low = calc_contraction(glo, B_up);
  return calc_contraction(B_low, B_up); // Bsq
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DEntropy::get_BiSi_Exact(const vec<CCTK_REAL, 3> &Bvec,
                              const vec<CCTK_REAL, 3> &mom) const {

  return Bvec(X) * mom(X) + Bvec(Y) * mom(Y) + Bvec(Z) * mom(Z); // BiSi
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 3>
c2p_1DEntropy::get_WLorentz_vsq_bsq_Seeds(const vec<CCTK_REAL, 3> &B_up,
                                          const vec<CCTK_REAL, 3> &v_up,
                                          const smat<CCTK_REAL, 3> &glo) const {
  vec<CCTK_REAL, 3> v_low = calc_contraction(glo, v_up);
  CCTK_REAL vsq = calc_contraction(v_low, v_up);
  CCTK_REAL VdotB = calc_contraction(v_low, B_up);
  CCTK_REAL VdotBsq = VdotB * VdotB;
  CCTK_REAL Bsq = get_Bsq_Exact(B_up, glo);

  CCTK_REAL w_lor = 1. / sqrt(1. - vsq);
  CCTK_REAL bsq = ((Bsq) / (w_lor * w_lor)) + VdotBsq;
  vec<CCTK_REAL, 3> w_vsq_bsq{w_lor, vsq, bsq};

  return w_vsq_bsq; //{w_lor, vsq, bsq}
}

/* Called by 1DEntropy */

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_1DEntropy::xEntropyToPrim(CCTK_REAL xEntropy_Sol, CCTK_REAL Ssq,
                              CCTK_REAL Bsq, CCTK_REAL BiSi,
                              const EOSType *eos_3p, prim_vars &pv,
                              const cons_vars &cv,
                              const smat<CCTK_REAL, 3> &gup,
                              const smat<CCTK_REAL, 3> &glo) const {
  // Density, entropy, Ye
  pv.rho = xEntropy_Sol;
  pv.entropy = cv.DEnt / cv.dens;
  pv.Ye = cv.DYe / cv.dens;

  // Lorentz factor
  pv.w_lor = cv.dens / xEntropy_Sol;

  // Pressure and epsilon
  pv.press =
      eos_3p->press_from_valid_rho_kappa_ye(xEntropy_Sol, pv.entropy, pv.Ye);
  pv.eps = eos_3p->eps_from_valid_rho_kappa_ye(xEntropy_Sol, pv.entropy, pv.Ye);

  // Temperature
  pv.temperature =
      eos_3p->temp_from_valid_rho_eps_ye(xEntropy_Sol, pv.eps, pv.Ye);

  // Taken from WZ2Prim (2DNRNoble)
  // Z_Sol = rho * h * w_lor * w_lor
  CCTK_REAL Z_Sol =
      (xEntropy_Sol + pv.eps * xEntropy_Sol + pv.press) * pv.w_lor * pv.w_lor;

  // TODO: Debug code to capture v>1,
  // remove soon
  if (use_zprim) {

    CCTK_REAL zx = pv.w_lor *
                   (gup(X, X) * cv.mom(X) + gup(X, Y) * cv.mom(Y) +
                    gup(X, Z) * cv.mom(Z)) /
                   (Z_Sol + Bsq);
    zx += pv.w_lor * BiSi * cv.dBvec(X) / (Z_Sol * (Z_Sol + Bsq));

    CCTK_REAL zy = pv.w_lor *
                   (gup(X, Y) * cv.mom(X) + gup(Y, Y) * cv.mom(Y) +
                    gup(Y, Z) * cv.mom(Z)) /
                   (Z_Sol + Bsq);
    zy += pv.w_lor * BiSi * cv.dBvec(Y) / (Z_Sol * (Z_Sol + Bsq));

    CCTK_REAL zz = pv.w_lor *
                   (gup(X, Z) * cv.mom(X) + gup(Y, Z) * cv.mom(Y) +
                    gup(Z, Z) * cv.mom(Z)) /
                   (Z_Sol + Bsq);
    zz += pv.w_lor * BiSi * cv.dBvec(Z) / (Z_Sol * (Z_Sol + Bsq));

    CCTK_REAL zx_down = glo(X, X) * zx + glo(X, Y) * zy + glo(X, Z) * zz;
    CCTK_REAL zy_down = glo(X, Y) * zx + glo(Y, Y) * zy + glo(Y, Z) * zz;
    CCTK_REAL zz_down = glo(X, Z) * zx + glo(Y, Z) * zy + glo(Z, Z) * zz;

    CCTK_REAL Zsq = zx * zx_down + zy * zy_down + zz * zz_down;

    CCTK_REAL SafeLor = sqrt(1.0 + Zsq);

    pv.vel(X) = zx / SafeLor;
    pv.vel(Y) = zy / SafeLor;
    pv.vel(Z) = zz / SafeLor;

    pv.w_lor = SafeLor;

  } else {

    pv.vel(X) = (gup(X, X) * cv.mom(X) + gup(X, Y) * cv.mom(Y) +
                 gup(X, Z) * cv.mom(Z)) /
                (Z_Sol + Bsq);
    pv.vel(X) += BiSi * cv.dBvec(X) / (Z_Sol * (Z_Sol + Bsq));

    pv.vel(Y) = (gup(X, Y) * cv.mom(X) + gup(Y, Y) * cv.mom(Y) +
                 gup(Y, Z) * cv.mom(Z)) /
                (Z_Sol + Bsq);
    pv.vel(Y) += BiSi * cv.dBvec(Y) / (Z_Sol * (Z_Sol + Bsq));

    pv.vel(Z) = (gup(X, Z) * cv.mom(X) + gup(Y, Z) * cv.mom(Y) +
                 gup(Z, Z) * cv.mom(Z)) /
                (Z_Sol + Bsq);
    pv.vel(Z) += BiSi * cv.dBvec(Z) / (Z_Sol * (Z_Sol + Bsq));
  }

  pv.Bvec = cv.dBvec;

  const vec<CCTK_REAL, 3> Elow = calc_cross_product(pv.Bvec, pv.vel);
  pv.E = calc_contraction(gup, Elow);
}

// See Appendix A.4 of https://arxiv.org/pdf/1112.0568
template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DEntropy::funcRoot_1DEntropy(CCTK_REAL Ssq, CCTK_REAL Bsq, CCTK_REAL BiSi,
                                  CCTK_REAL x, const EOSType *eos_3p,
                                  const cons_vars &cv) const {

  // We already divided dens, DEnt and
  // DYe by sqrt(gamma)

  const CCTK_REAL ent_loc = cv.DEnt / cv.dens;
  const CCTK_REAL ye_loc = cv.DYe / cv.dens;

  // Compute h using entropy
  const CCTK_REAL press_loc =
      eos_3p->press_from_valid_rho_kappa_ye(x, ent_loc, ye_loc);

  const CCTK_REAL eps_loc =
      eos_3p->eps_from_valid_rho_kappa_ye(x, ent_loc, ye_loc);

  // Compute (A60) using
  // W = rho*h*lorentz*lorentz
  const CCTK_REAL lor_loc = cv.dens / x;
  const CCTK_REAL h_loc = (1.0 + eps_loc + press_loc / x);
  const CCTK_REAL W = x * h_loc * lor_loc * lor_loc;

  // Compute (A61)
  const CCTK_REAL Sfsq =
      (W * W * Ssq + BiSi * BiSi * (Bsq + 2.0 * W)) / ((W + Bsq) * (W + Bsq));

  // Compute (A62)
  const CCTK_REAL rho =
      cv.dens / sqrt(1.0 + Sfsq / (cv.dens * cv.dens * h_loc * h_loc));

  return x - rho;
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_1DEntropy::solve(const EOSType *eos_3p, prim_vars &pv, cons_vars &cv,
                     const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const {

  // ROOTSTAT status = ROOTSTAT::SUCCESS;
  rep.iters = 0;
  rep.adjust_cons = false;
  rep.set_atmo = false;
  rep.status = c2p_report::SUCCESS;

  /* Check validity of the 3-metric and compute its inverse */
  const CCTK_REAL spatial_detg = calc_det(glo);
  const CCTK_REAL sqrt_detg = sqrt(spatial_detg);

  // if ((!isfinite(sqrt_detg)) || (sqrt_detg <= 0)) {
  //  rep.set_invalid_detg(sqrt_detg);
  //  set_to_nan(pv, cv);
  //  return;
  //}

  // Check positive definiteness of spatial metric
  // Sylvester's criterion, see
  // https://en.wikipedia.org/wiki/Sylvester%27s_criterion
  const bool minor1{glo(X, X) > 0.0};
  const bool minor2{glo(X, X) * glo(Y, Y) - glo(X, Y) * glo(X, Y) > 0.0};
  const bool minor3{spatial_detg > 0.0};

  if (!(minor1 && minor2 && minor3)) {
    rep.set_invalid_detg(sqrt_detg);
    set_to_nan(pv, cv);
    return;
  }

  const smat<CCTK_REAL, 3> gup = calc_inv(glo, spatial_detg);

  /* Undensitize the conserved vars */
  /* Make sure to return densitized values later on! */
  cv.dens /= sqrt_detg;
  cv.tau /= sqrt_detg;
  cv.mom /= sqrt_detg;
  cv.dBvec /= sqrt_detg;
  cv.DYe /= sqrt_detg;
  cv.DEnt /= sqrt_detg;

  // if (cv.dens <= atmo.rho_cut) {
  //  rep.set_atmo_set();
  //  pv.Bvec = cv.dBvec;
  //  atmo.set(pv, cv, glo);
  //  return;
  //}

  const CCTK_REAL Ssq = get_Ssq_Exact(cv.mom, gup);
  const CCTK_REAL Bsq = get_Bsq_Exact(cv.dBvec, glo);
  const CCTK_REAL BiSi = get_BiSi_Exact(cv.dBvec, cv.mom);

  // TODO: Is this check really necessary?
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

  // See Appendix A.4 of https://arxiv.org/pdf/1112.0568
  // Compute (A59)
  CCTK_REAL a = cv.dens / sqrt(1.0 + Ssq / (cv.dens * cv.dens));
  CCTK_REAL b = cv.dens;
  auto fn = [&](auto x) {
    return funcRoot_1DEntropy(Ssq, Bsq, BiSi, x, eos_3p, cv);
  };

  // Important!
  // Algo::brent terminates if the following accuracy is achieved
  // abs(x - y) <= eps * min(abs(x), abs(y)),
  // where x and y are the values of the bracket and
  // eps = std::ldexp(1, -minbits) = 1 * 2^{-minbits}
  // This should probably be changed in Algo::brent

  // We want to set the tolerance to its correct parameter
  const CCTK_REAL log2 = std::log(2.0);
  const CCTK_INT minbits = int(abs(std::log(tolerance)) / log2);
  const CCTK_REAL tolerance_0 = std::ldexp(double(1.0), -minbits);
  // Old code:
  // const CCTK_INT minbits = std::numeric_limits<CCTK_REAL>::digits - 4;
  const CCTK_INT maxiters = maxIterations;

  auto result = Algo::brent(fn, a, b, minbits, maxiters, rep.iters);

  CCTK_REAL xEntropy_Sol = 0.5 * (result.first + result.second);

  xEntropyToPrim(xEntropy_Sol, Ssq, Bsq, BiSi, eos_3p, pv, cv, gup, glo);

  // Check solution and calculate primitives
  //  if (rep.iters < maxiters && abs(fn(xEntropy_Sol)) < tolerance) {
  /*
  if (abs(result.first - result.second) <=
      tolerance_0 * min(abs(result.first), abs(result.second))) {
    rep.status = c2p_report::SUCCESS;
    status = ROOTSTAT::SUCCESS;
  } else {
    // set status to root not converged
    rep.set_root_conv();
    status = ROOTSTAT::NOT_CONVERGED;
  }
  */

  // Error out if rho is negative or zero
  if (pv.rho <= 0.0) {
    // set status to rho is out of range
    rep.set_range_rho(cv.dens, pv.rho);
    cv.dens *= sqrt_detg;
    cv.tau *= sqrt_detg;
    cv.mom *= sqrt_detg;
    cv.dBvec *= sqrt_detg;
    cv.DYe *= sqrt_detg;
    cv.DEnt *= sqrt_detg;
    return;
  }

  if (abs(result.first - result.second) >
      tolerance_0 * min(abs(result.first), abs(result.second))) {

    // check primitives against conservatives
    cons_vars cv_check;
    cv_check.from_prim(pv, glo);
    /* Undensitize the conserved vars */
    cv_check.dens /= sqrt_detg;
    cv_check.tau /= sqrt_detg;
    cv_check.mom /= sqrt_detg;
    cv_check.dBvec /= sqrt_detg;
    cv_check.DYe /= sqrt_detg;
    cv_check.DEnt /= sqrt_detg;

    CCTK_REAL small = 1e-50;

    // note that we don't compute the error in
    // tau as we make use of the conserved entropy
    // for the inversion
    CCTK_REAL max_error = sqrt(
        max({pow((cv_check.dens - cv.dens) / (cv.dens + small), 2.0),
             pow((cv_check.mom(0) - cv.mom(0)) / (cv.mom(0) + small), 2.0),
             pow((cv_check.mom(1) - cv.mom(1)) / (cv.mom(1) + small), 2.0),
             pow((cv_check.mom(2) - cv.mom(2)) / (cv.mom(2) + small), 2.0)}));

    // reject only if mismatch in conservatives is inappropriate, else accept
    if (max_error >= cons_error) {
      // set status to root not converged
      rep.set_root_conv();
      // status = ROOTSTAT::NOT_CONVERGED;
      cv.dens *= sqrt_detg;
      cv.tau *= sqrt_detg;
      cv.mom *= sqrt_detg;
      cv.dBvec *= sqrt_detg;
      cv.DYe *= sqrt_detg;
      cv.DEnt *= sqrt_detg;
      return;
    }
  }

  // Lower velocity
  const vec<CCTK_REAL, 3> v_low = calc_contraction(glo, pv.vel);
  // Computing b^t : this is b^0 * alp
  const CCTK_REAL bst = pv.w_lor * calc_contraction(pv.Bvec, v_low);
  // Computing b^mu b_mu
  const CCTK_REAL bs2 = (Bsq + bst * bst) / (pv.w_lor * pv.w_lor);
  // Recompute tau
  cv.tau = (pv.w_lor * pv.w_lor * (pv.rho * (1.0 + pv.eps) + pv.press + bs2) -
            (pv.press + 0.5 * bs2) - bst * bst) -
           cv.dens;

  // set to atmo if computed rho is below floor density
  if (pv.rho < atmo.rho_cut) {
    rep.set_atmo_set();
    atmo.set(pv, cv, glo);
    return;
  }

  c2p::prims_floors_and_ceilings(eos_3p, pv, cv, glo, rep);

  // Recompute cons if prims have been adjusted
  if (rep.adjust_cons) {
    cv.from_prim(pv, glo);
  } else {
    /* Densitize the conserved vars again*/
    cv.dens *= sqrt_detg;
    cv.tau *= sqrt_detg;
    cv.mom *= sqrt_detg;
    cv.dBvec *= sqrt_detg;
    cv.DYe *= sqrt_detg;
    cv.DEnt *= sqrt_detg;
  }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_1DEntropy::set_to_nan(prim_vars &pv, cons_vars &cv) const {
  pv.set_to_nan();
  cv.set_to_nan();
}

/* Destructor */
CCTK_HOST CCTK_DEVICE
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DEntropy::~c2p_1DEntropy() {
  // How to destruct properly a vector?
}
} // namespace Con2PrimFactory

#endif
