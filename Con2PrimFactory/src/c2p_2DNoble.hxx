#ifndef C2P_2DNOBLE_HXX
#define C2P_2DNOBLE_HXX

#include "c2p.hxx"

namespace Con2PrimFactory {

using namespace std;

class c2p_2DNoble : public c2p {
public:
  /* Some attributes */
  CCTK_REAL GammaIdealFluid;

  /* Constructor */
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_2DNoble(
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

  /* Called by 2DNoble */
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Z_Seed(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL press,
             CCTK_REAL w_lor) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Press_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq, cons_vars &cv) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq, cons_vars &cv) const;
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  WZ2Prim(CCTK_REAL Z_Sol, CCTK_REAL vsq_Sol, CCTK_REAL Bsq, CCTK_REAL BiSi,
          EOSType &eos_th, prim_vars &pv, cons_vars cv,
          const smat<CCTK_REAL, 3> &gup) const;
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  solve(EOSType &eos_th, prim_vars &pv, prim_vars &pv_seeds, cons_vars cv,
        const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const;

  /* Destructor */
  CCTK_HOST CCTK_DEVICE ~c2p_2DNoble();
};

/* Constructor */
template <typename EOSType>
CCTK_HOST
    CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_2DNoble::c2p_2DNoble(
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
c2p_2DNoble::get_Ssq_Exact(vec<CCTK_REAL, 3> &mom_low,
                           const smat<CCTK_REAL, 3> &gup) const {
  vec<CCTK_REAL, 3> mom_up = calc_contraction(gup, mom_low);
  CCTK_REAL Ssq = calc_contraction(mom_low, mom_up);

  if ((Ssq < 0.) && (fabs(Ssq) < 1.0e-13)) {
    Ssq = fabs(Ssq);
  }
  return Ssq;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_Bsq_Exact(vec<CCTK_REAL, 3> &B_up,
                           const smat<CCTK_REAL, 3> &glo) const {
  vec<CCTK_REAL, 3> B_low = calc_contraction(glo, B_up);
  return calc_contraction(B_low, B_up); // Bsq
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_BiSi_Exact(vec<CCTK_REAL, 3> &Bvec,
                            vec<CCTK_REAL, 3> &mom) const {
  return calc_contraction(mom, Bvec); // BiS^i
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 3>
c2p_2DNoble::get_WLorentz_vsq_bsq_Seeds(vec<CCTK_REAL, 3> &B_up,
                                        vec<CCTK_REAL, 3> &v_up,
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

/* Called by 2dNRNoble */

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_Z_Seed(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL press,
                        CCTK_REAL w_lor) const {
  return (rho + eps * rho + press) * w_lor * w_lor;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_Press_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq,
                                cons_vars &cv) const {
  return ((Z * (1.0 - Vsq) - cv.dens * sqrt(1.0 - Vsq)) *
          (GammaIdealFluid - 1.0) / (GammaIdealFluid));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq) const {
  return ((1.0 - Vsq) * (GammaIdealFluid - 1.0) / GammaIdealFluid);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq,
                                 cons_vars &cv) const {
  return ((-Z + cv.dens / (2.0 * sqrt(1.0 - Vsq))) * (GammaIdealFluid - 1.0) /
          GammaIdealFluid);
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_2DNoble::WZ2Prim(CCTK_REAL Z_Sol, CCTK_REAL vsq_Sol, CCTK_REAL Bsq,
                     CCTK_REAL BiSi, EOSType &eos_th, prim_vars &pv,
                     cons_vars cv, const smat<CCTK_REAL, 3> &gup) const {
  CCTK_REAL W_Sol = 1.0 / sqrt(1.0 - vsq_Sol);

  pv.rho = cv.dens / W_Sol;

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

  pv.w_lor = W_Sol;

  pv.eps = (Z_Sol * (1. - vsq_Sol) / pv.rho - 1.0) / GammaIdealFluid;

  pv.Ye = cv.dYe / cv.dens;

  pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

  pv.Bvec = cv.dBvec;

  const vec<CCTK_REAL, 3> Elow = calc_cross_product(pv.Bvec, pv.vel);
  pv.E = calc_contraction(gup, Elow);
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_2DNoble::solve(EOSType &eos_th, prim_vars &pv, prim_vars &pv_seeds,
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
  CCTK_REAL vsq_seed = w_vsq_bsq(1);
  // const CCTK_REAL bsq = w_vsq_bsq(2); // small b^2

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

  CCTK_REAL eps_min =
      eos_th.range_eps_from_valid_rho_ye(pv_seeds.rho, pv_seeds.Ye).min;
  CCTK_REAL eps_last = max({pv_seeds.eps, eps_min});

  /* get pressure seed from updated pv_seeds.rho */
  pv_seeds.press =
      eos_th.press_from_valid_rho_eps_ye(pv_seeds.rho, eps_last, pv_seeds.Ye);

  /* get Z seed */
  CCTK_REAL Z_Seed =
      get_Z_Seed(pv_seeds.rho, eps_last, pv_seeds.press, pv_seeds.w_lor);

  /* initialize unknowns for c2p, Z and vsq: */
  CCTK_REAL x[2];
  CCTK_REAL x_old[2];
  x[0] = fabs(Z_Seed);
  x[1] = vsq_seed;

  //      printf("x[0], x[1]: %f, %f \n", x[0], x[1]);

  /* initialize old values */
  x_old[0] = x[0];
  x_old[1] = x[1];

  /* Start Recovery with 2D NR Solver */
  const CCTK_INT n = 2;
  const CCTK_REAL dv = (1. - 1.e-15);
  CCTK_REAL dx[n];
  CCTK_REAL fjac[n][n];
  CCTK_REAL resid[n];

  CCTK_REAL errx = 1.;
  CCTK_REAL df = 1.;
  CCTK_REAL f = 1.;

  CCTK_INT k;
  for (k = 1; k <= maxIterations; k++) {

    /* Expressions for the jacobian are adapted from the Noble C2P
    implementation in the Spritz code. As the analytical form of the equations
    is known, the Newton-Raphson step can be computed explicitly */

    const CCTK_REAL Z = x[0];
    const CCTK_REAL invZ = 1.0 / Z;
    const CCTK_REAL Vsq = x[1];
    const CCTK_REAL Sdotn = -(cv.tau + cv.dens);
    const CCTK_REAL p_tmp = get_Press_funcZVsq(Z, Vsq, cv);
    const CCTK_REAL dPdvsq = get_dPdVsq_funcZVsq(Z, Vsq, cv);
    const CCTK_REAL dPdZ = get_dPdZ_funcZVsq(Z, Vsq);

    fjac[0][0] = -2 * (Vsq + BiSi * BiSi * invZ * invZ * invZ) * (Bsq + Z);
    fjac[0][1] = -(Bsq + Z) * (Bsq + Z);
    fjac[1][0] = -1.0 + dPdZ - BiSi * BiSi * invZ * invZ * invZ;
    fjac[1][1] = -0.5 * Bsq + dPdvsq;

    resid[0] = Ssq - Vsq * (Bsq + Z) * (Bsq + Z) +
               BiSi * BiSi * invZ * invZ * (Bsq + Z + Z);
    resid[1] = -Sdotn - 0.5 * Bsq * (1.0 + Vsq) +
               0.5 * BiSi * BiSi * invZ * invZ - Z + p_tmp;

    const CCTK_REAL detjac =
        (Bsq + Z) *
        (fjac[1][0] * (Bsq + Z) +
         (Bsq - 2.0 * dPdvsq) * (BiSi * BiSi * invZ * invZ + Vsq * Z) * invZ);
    const CCTK_REAL detjac_inv = 1.0 / detjac;

    dx[0] = -(fjac[1][1] * resid[0] - fjac[0][1] * resid[1]) * detjac_inv;
    dx[1] = -(-fjac[1][0] * resid[0] + fjac[0][0] * resid[1]) * detjac_inv;

    df = -resid[0] * resid[0] - resid[1] * resid[1];
    f = -0.5 * (df);

    /* save old values before calculating the new */
    errx = 0.;
    x_old[0] = x[0];
    x_old[1] = x[1];

    // make the newton step
    x[0] += dx[0];
    x[1] += dx[1];

    // calculate the convergence criterion
    errx = (x[0] == 0.) ? fabs(dx[0]) : fabs(dx[0] / x[0]);

    /* make sure that the new x[] is physical */
    if (x[0] < 0.0) {
      x[0] = fabs(x[0]);
    } else {
      if (x[0] > 1e20) {
        x[0] = x_old[0];
      }
    }

    if (x[1] < 0.0) {
      x[1] = 0.0;
    }

    else {
      if (x[1] >= 1.0) {
        x[1] = dv;
      }
    }
  }

  // storing number of iterations taken to find the root
  rep.iters = k;

  if (fabs(errx) <= tolerance) {
    rep.status = c2p_report::SUCCESS;
    status = ROOTSTAT::SUCCESS;
  } else {
    // set status to root not converged
    rep.set_root_conv();
    status = ROOTSTAT::NOT_CONVERGED;
  }

  // Check for bad untrapped divergences
  if ((!isfinite(f)) || (!isfinite(df))) {
    rep.set_root_bracket();
    // TODO: currently, divergences in f and df are marked as failed bracketing.
    // Change this.
  }

  /* Calculate primitives from Z and W */
  CCTK_REAL Z_Sol = x[0];
  CCTK_REAL vsq_Sol = x[1];

  /* Write prims if C2P succeeded */
  WZ2Prim(Z_Sol, vsq_Sol, Bsq, BiSi, eos_th, pv, cv, gup);

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
    rep.set_range_eps(rgeps.min); // sets adjust_cons to true by default
  }

  // TODO: check validity for Ye

  // check if computed velocities are within the specified limit
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
c2p_2DNoble::set_to_nan(prim_vars &pv, cons_vars &cv) const {
  pv.set_to_nan();
  cv.set_to_nan();
}

/* Destructor */
CCTK_HOST CCTK_DEVICE
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_2DNoble::~c2p_2DNoble() {
  // How to destruct properly a vector?
}
} // namespace Con2PrimFactory

#endif
