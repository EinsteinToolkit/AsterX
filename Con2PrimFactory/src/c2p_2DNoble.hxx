#ifndef C2P_2DNOBLE_HXX
#define C2P_2DNOBLE_HXX

#include "c2p.hxx"

namespace Con2PrimFactory {

using namespace std;

class c2p_2DNoble : public c2p {
public:
  /* Some attributes */
  CCTK_REAL GammaIdealFluid;
  CCTK_REAL Zmin;

  /* Constructor */
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_2DNoble(
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
  // TODO: Debug function to capture v>1,
  // remove soon
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 3>
  getZ_WLorentz_vsq_bsq_Seeds(const vec<CCTK_REAL, 3> &B_up,
                              const vec<CCTK_REAL, 3> &z_up,
                              const smat<CCTK_REAL, 3> &glo) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_to_nan(prim_vars &pv, cons_vars &cv) const;

  /* Called by 2DNoble */
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Z_Seed(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL press,
             CCTK_REAL w_lor) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Press_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq, const cons_vars &cv) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq, const cons_vars &cv) const;
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  WZ2Prim(CCTK_REAL Z_Sol, CCTK_REAL vsq_Sol, CCTK_REAL Bsq, CCTK_REAL BiSi,
          const EOSType *eos_3p, prim_vars &pv, const cons_vars &cv,
          const smat<CCTK_REAL, 3> &gup, const smat<CCTK_REAL, 3> &glo) const;
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  solve(const EOSType *eos_3p, prim_vars &pv, prim_vars &pv_seeds,
        cons_vars &cv, const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const;

  /* Destructor */
  CCTK_HOST CCTK_DEVICE ~c2p_2DNoble();
};

/* Constructor */
template <typename EOSType>
CCTK_HOST
    CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_2DNoble::c2p_2DNoble(
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
  Zmin = atmo.rho_atmo;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_Ssq_Exact(const vec<CCTK_REAL, 3> &mom_low,
                           const smat<CCTK_REAL, 3> &gup) const {
  vec<CCTK_REAL, 3> mom_up = calc_contraction(gup, mom_low);
  CCTK_REAL Ssq = calc_contraction(mom_low, mom_up);

  return Ssq;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_Bsq_Exact(const vec<CCTK_REAL, 3> &B_up,
                           const smat<CCTK_REAL, 3> &glo) const {
  vec<CCTK_REAL, 3> B_low = calc_contraction(glo, B_up);
  return calc_contraction(B_low, B_up); // Bsq
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_BiSi_Exact(const vec<CCTK_REAL, 3> &Bvec,
                            const vec<CCTK_REAL, 3> &mom) const {
  return calc_contraction(mom, Bvec); // BiS^i
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 3>
c2p_2DNoble::get_WLorentz_vsq_bsq_Seeds(const vec<CCTK_REAL, 3> &B_up,
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

// TODO: Debug function to capture v>1,
// remove soon
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 3>
c2p_2DNoble::getZ_WLorentz_vsq_bsq_Seeds(const vec<CCTK_REAL, 3> &B_up,
                                         const vec<CCTK_REAL, 3> &z_up,
                                         const smat<CCTK_REAL, 3> &glo) const {
  vec<CCTK_REAL, 3> z_low = calc_contraction(glo, z_up);
  CCTK_REAL zsq = calc_contraction(z_low, z_up);

  CCTK_REAL w_lor = sqrt(1. + zsq);
  CCTK_REAL vsq = min(zsq / w_lor / w_lor, 1. - 1.e-15);

  CCTK_REAL VdotB = calc_contraction(z_low, B_up) / w_lor;
  CCTK_REAL VdotBsq = VdotB * VdotB;
  CCTK_REAL Bsq = get_Bsq_Exact(B_up, glo);

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
                                const cons_vars &cv) const {
  return ((Z * (1.0 - Vsq) - cv.dens * sqrt(1.0 - Vsq)) *
          (GammaIdealFluid - 1.0) / (GammaIdealFluid));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq) const {
  return ((1.0 - Vsq) * (GammaIdealFluid - 1.0) / GammaIdealFluid);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq,
                                 const cons_vars &cv) const {
  return ((-Z + cv.dens / (2.0 * sqrt(1.0 - Vsq))) * (GammaIdealFluid - 1.0) /
          GammaIdealFluid);
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_2DNoble::WZ2Prim(CCTK_REAL Z_Sol, CCTK_REAL vsq_Sol, CCTK_REAL Bsq,
                     CCTK_REAL BiSi, const EOSType *eos_3p, prim_vars &pv,
                     const cons_vars &cv, const smat<CCTK_REAL, 3> &gup,
                     const smat<CCTK_REAL, 3> &glo) const {
  CCTK_REAL W_Sol = 1.0 / sqrt(1.0 - vsq_Sol);

  pv.rho = cv.dens / W_Sol;

  // TODO: Debug code to capture v>1,
  // remove soon
  if (use_zprim) {

    CCTK_REAL zx = W_Sol *
                   (gup(X, X) * cv.mom(X) + gup(X, Y) * cv.mom(Y) +
                    gup(X, Z) * cv.mom(Z)) /
                   (Z_Sol + Bsq);
    zx += W_Sol * BiSi * cv.dBvec(X) / (Z_Sol * (Z_Sol + Bsq));

    CCTK_REAL zy = W_Sol *
                   (gup(X, Y) * cv.mom(X) + gup(Y, Y) * cv.mom(Y) +
                    gup(Y, Z) * cv.mom(Z)) /
                   (Z_Sol + Bsq);
    zy += W_Sol * BiSi * cv.dBvec(Y) / (Z_Sol * (Z_Sol + Bsq));

    CCTK_REAL zz = W_Sol *
                   (gup(X, Z) * cv.mom(X) + gup(Y, Z) * cv.mom(Y) +
                    gup(Z, Z) * cv.mom(Z)) /
                   (Z_Sol + Bsq);
    zz += W_Sol * BiSi * cv.dBvec(Z) / (Z_Sol * (Z_Sol + Bsq));

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

    pv.w_lor = W_Sol;
  }

  // pv.eps = (Z_Sol * (1. - vsq_Sol) / pv.rho - 1.0) / GammaIdealFluid;
  pv.eps = (Z_Sol / pv.w_lor / pv.w_lor / pv.rho - 1.0) / GammaIdealFluid;

  pv.Ye = cv.DYe / cv.dens;

  pv.press = eos_3p->press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

  pv.temperature = eos_3p->temp_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

  pv.entropy = eos_3p->kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

  pv.Bvec = cv.dBvec;

  const vec<CCTK_REAL, 3> Elow = calc_cross_product(pv.Bvec, pv.vel);
  pv.E = calc_contraction(gup, Elow);
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_2DNoble::solve(const EOSType *eos_3p, prim_vars &pv, prim_vars &pv_seeds,
                   cons_vars &cv, const smat<CCTK_REAL, 3> &glo,
                   c2p_report &rep) const {

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

  // compute primitive B seed from conserved B of current time step for better
  // guess
  pv_seeds.Bvec = cv.dBvec;

  const CCTK_REAL Ssq = get_Ssq_Exact(cv.mom, gup);
  const CCTK_REAL Bsq = get_Bsq_Exact(pv_seeds.Bvec, glo);
  const CCTK_REAL BiSi = get_BiSi_Exact(pv_seeds.Bvec, cv.mom);

  vec<CCTK_REAL, 3> w_vsq_bsq;
  CCTK_REAL vsq_seed;

  // TODO: Debug code to capture v>1,
  // remove soon
  if (use_zprim) {

    vec<CCTK_REAL, 3> zvec = pv_seeds.vel * pv_seeds.w_lor;

    w_vsq_bsq = getZ_WLorentz_vsq_bsq_Seeds(
        pv_seeds.Bvec, zvec, glo); // this also recomputes pv_seeds.w_lor

    pv_seeds.w_lor = w_vsq_bsq(0);
    vsq_seed = w_vsq_bsq(1);

  } else {

    w_vsq_bsq =
        get_WLorentz_vsq_bsq_Seeds(pv_seeds.Bvec, pv_seeds.vel,
                                   glo); // this also recomputes pv_seeds.w_lor
    pv_seeds.w_lor = w_vsq_bsq(0);
    vsq_seed = w_vsq_bsq(1);
  }

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

  /* update rho seed from cv and wlor */
  // rho consistent with cv.rho should be better guess than rho from last
  // timestep
  pv_seeds.rho = cv.dens / pv_seeds.w_lor;

  CCTK_REAL eps_last = max({pv_seeds.eps, atmo.eps_atmo});

  /* get pressure seed from updated pv_seeds.rho */
  pv_seeds.press =
      eos_3p->press_from_valid_rho_eps_ye(pv_seeds.rho, eps_last, pv_seeds.Ye);

  /* get Z seed */
  CCTK_REAL Z_Seed =
      get_Z_Seed(pv_seeds.rho, eps_last, pv_seeds.press, pv_seeds.w_lor);

  /* initialize unknowns for c2p, Z and vsq: */
  CCTK_REAL x[2];
  CCTK_REAL x_old[2];
  x[0] = fabs(Z_Seed);
  x[1] = vsq_seed;

  /* initialize old values */
  x_old[0] = x[0];
  x_old[1] = x[1];

  /* Start Recovery with 2D NR Solver */
  constexpr CCTK_INT n = 2;
  constexpr CCTK_REAL dv = (1. - 1.e-10);
  constexpr CCTK_REAL dw = 1. / (1. - dv);

  CCTK_REAL dx[n];
  CCTK_REAL fjac[n][n];
  CCTK_REAL resid[n];

  CCTK_REAL errx = 1.;
  CCTK_REAL df = 1.;
  CCTK_REAL f = 1.;

  /* make sure that x[] is physical */
  if (x[1] < 0.0) {
    x[1] = 0.0;
  }

  else {
    if (x[1] >= dv) {
      x[0] *= dw * (1.0 - x[1]);
      x[1] = dv;
    }
  }

  if (x[0] < Zmin) {
    x[0] = Zmin;
  } else {
    if (x[0] > 1e20) {
      x[0] = x_old[0];
    }
  }

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

    /* make sure that the new x[] is physical */
    if (x[1] < 0.0) {
      x[1] = 0.0;
    }

    else {
      if (x[1] >= dv) {
        x[0] *= dw * (1.0 - x[1]);
        x[1] = dv;
      }
    }

    if (x[0] < Zmin) {
      x[0] = Zmin;
    } else {
      if (x[0] > 1e20) {
        x[0] = x_old[0];
      }
    }

    // calculate the convergence criterion
    errx = (x[0] == 0.) ? fabs(dx[0]) : fabs(dx[0] / x[0]);

    if (fabs(errx) <= tolerance) {
      break;
    }
  }

  // storing number of iterations taken to find the root
  rep.iters = k;

  // if (fabs(errx) <= tolerance) {
  // rep.status = c2p_report::SUCCESS;
  // status = ROOTSTAT::SUCCESS;
  //} else {
  // set status to root not converged
  // rep.set_root_conv();
  // status = ROOTSTAT::NOT_CONVERGED;
  //}

  if (fabs(errx) > tolerance) {
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
  WZ2Prim(Z_Sol, vsq_Sol, Bsq, BiSi, eos_3p, pv, cv, gup, glo);

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

  // Error out if eps is negative or zero
  if (pv.eps <= 0.0) {
    // set status to eps is out of range
    rep.set_range_eps(pv.eps);
    cv.dens *= sqrt_detg;
    cv.tau *= sqrt_detg;
    cv.mom *= sqrt_detg;
    cv.dBvec *= sqrt_detg;
    cv.DYe *= sqrt_detg;
    cv.DEnt *= sqrt_detg;
    return;
  }

  // Conserved entropy must be consistent with new prims
  cv.DEnt = cv.dens * pv.entropy;

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
    /* If not adjusted, densitize back the original conserved vars */
    cv.dens *= sqrt_detg;
    cv.tau *= sqrt_detg;
    cv.mom *= sqrt_detg;
    cv.dBvec *= sqrt_detg;
    cv.DYe *= sqrt_detg;
    cv.DEnt *= sqrt_detg;
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
