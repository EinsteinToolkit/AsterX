#ifndef C2P_2DNOBLE_HXX
#define C2P_2DNOBLE_HXX

#include "c2p.hxx"

namespace Con2PrimFactory {

using namespace std;

class c2p_2DNoble : public c2p {
public:
  /* Some attributes */
  CCTK_REAL Zmin;

  /* Constructor */
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_2DNoble(
      const EOSType *eos_3p, const atmosphere &atm, CCTK_INT maxIter,
      CCTK_REAL tol, CCTK_REAL alp_thresh_in, CCTK_REAL vwlim, CCTK_REAL B_lim,
      CCTK_REAL rho_BH_in, CCTK_REAL eps_BH_in, CCTK_REAL vwlim_BH_in,
      CCTK_REAL sigma_max_in, CCTK_REAL inv_beta_max_in, bool ye_len,
      bool use_z, bool use_temperature, bool use_pressure_atmo,
      bool soft_root_conv, CCTK_REAL soft_root_width_factor_in);

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
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE inline bool
  get_Press_funcZVsq(CCTK_REAL &press, CCTK_REAL &dPdZ,
                     CCTK_REAL &dPdVsq, CCTK_REAL Z, CCTK_REAL Vsq,
                     const EOSType *eos_3p, const cons_vars &cv,
                     CCTK_REAL *temp = nullptr) const;
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  WZ2Prim(CCTK_REAL Z_Sol, CCTK_REAL vsq_Sol, CCTK_REAL Bsq, CCTK_REAL BiSi,
          const EOSType *eos_3p, prim_vars &pv, CCTK_REAL &eps_raw,
          const cons_vars &cv,
          const smat<CCTK_REAL, 3> &gup, const smat<CCTK_REAL, 3> &glo) const;
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  solve(const EOSType *eos_3p, prim_vars &pv, prim_vars &pv_seeds,
        cons_vars &cv, const CCTK_REAL alp, const vec<CCTK_REAL, 3> &beta,
        const smat<CCTK_REAL, 3> &glo, c2p_report &rep,
        bool reject_nonpositive_eps = false) const;

  /* Destructor */
  CCTK_HOST CCTK_DEVICE ~c2p_2DNoble();
};

/* Constructor */
template <typename EOSType>
CCTK_HOST
    CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_2DNoble::c2p_2DNoble(
        const EOSType *eos_3p, const atmosphere &atm, CCTK_INT maxIter,
        CCTK_REAL tol, CCTK_REAL alp_thresh_in, CCTK_REAL vwlim,
        CCTK_REAL B_lim, CCTK_REAL rho_BH_in, CCTK_REAL eps_BH_in,
        CCTK_REAL vwlim_BH_in, CCTK_REAL sigma_max_in,
        CCTK_REAL inv_beta_max_in, bool ye_len, bool use_z,
        bool use_temperature, bool use_pressure_atmo, bool soft_root_conv,
        CCTK_REAL soft_root_width_factor_in) {

  // Base
  atmo = atm;
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
  sigma_max = sigma_max_in;
  inv_beta_max = inv_beta_max_in;
  ye_lenient = ye_len;
  use_zprim = use_z;
  use_temp = use_temperature;
  use_press_atmo = use_pressure_atmo;
  soft_root_convergence = soft_root_conv;
  soft_root_width_factor = fmax(CCTK_REAL(1.0), soft_root_width_factor_in);

  // Derived
  Zmin = eos_3p->rgrho.min;
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

template <typename EOSType>
CCTK_HOST CCTK_DEVICE inline bool
c2p_2DNoble::get_Press_funcZVsq(CCTK_REAL &press, CCTK_REAL &dPdZ,
                                CCTK_REAL &dPdVsq, CCTK_REAL Z,
                                CCTK_REAL Vsq, const EOSType *eos_3p,
                                const cons_vars &cv, CCTK_REAL *temp) const {
  if (!std::isfinite(Z) || Z <= 0.0 || !std::isfinite(Vsq) ||
      Vsq < 0.0 || Vsq >= 1.0 || !(cv.dens > 0.0))
    return false;
  const CCTK_REAL w_lor = 1.0 / sqrt(1.0 - Vsq);
  const CCTK_REAL rho = cv.dens / w_lor;
  const CCTK_REAL Ye =
      std::clamp(cv.DYe / cv.dens, eos_3p->rgye.min, eos_3p->rgye.max);
  const CCTK_REAL h = Z * (1.0 - Vsq) / rho;
  CCTK_REAL dpdrho, dpdeps;
  if (!eos_3p->press_derivs_from_rho_h_ye(press, dpdrho, dpdeps, rho, h, Ye, temp))
    return false;
  // Chain rule for Z = rho*h*W^2 and rho = D/W.
  const CCTK_REAL denom = 1.0 + dpdeps / rho;
  if (!std::isfinite(denom) || denom == 0.0)
    return false;
  dPdZ = (dpdeps / rho) * (1.0 - Vsq) / denom;
  dPdVsq = (-0.5 * cv.dens * w_lor * dpdrho -
            0.5 * dpdeps * (Z + press * w_lor * w_lor) / rho) / denom;
  return std::isfinite(press) && std::isfinite(dPdZ) && std::isfinite(dPdVsq);
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_2DNoble::WZ2Prim(CCTK_REAL Z_Sol, CCTK_REAL vsq_Sol, CCTK_REAL Bsq,
                     CCTK_REAL BiSi, const EOSType *eos_3p, prim_vars &pv,
                     CCTK_REAL &eps_raw, const cons_vars &cv,
                     const smat<CCTK_REAL, 3> &gup,
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

  pv.Ye = std::clamp(cv.DYe / cv.dens, eos_3p->rgye.min, eos_3p->rgye.max);
  if constexpr (std::is_same_v<EOSType, EOSX::eos_3p_idealgas>) {
    const CCTK_REAL h = Z_Sol / pv.w_lor / pv.w_lor / pv.rho;
    eos_3p->eps_from_rho_h_ye(pv.rho, h, pv.Ye, eps_raw);
  } else {
    const CCTK_REAL press_raw = Z_Sol + Bsq - cv.tau - cv.dens -
        0.5 * Bsq / (pv.w_lor * pv.w_lor) -
        0.5 * BiSi * BiSi / (Z_Sol * Z_Sol);
    eps_raw = (Z_Sol / (pv.w_lor * pv.w_lor) - press_raw) / pv.rho - 1.0;
  }
  const auto rgeps = eos_3p->range_eps_from_rho_ye(pv.rho, pv.Ye);
  pv.eps = std::clamp(eps_raw, rgeps.min, rgeps.max);

  pv.press = eos_3p->press_from_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

  pv.temperature = eos_3p->temp_from_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

  pv.entropy = eos_3p->kappa_from_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

  pv.Bvec = cv.dBvec;

  const vec<CCTK_REAL, 3> Elow = calc_cross_product(pv.Bvec, pv.vel);
  pv.E = calc_contraction(gup, Elow);
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_2DNoble::solve(const EOSType *eos_3p, prim_vars &pv, prim_vars &pv_seeds,
                   cons_vars &cv, const CCTK_REAL alp,
                   const vec<CCTK_REAL, 3> &beta, const smat<CCTK_REAL, 3> &glo,
                   c2p_report &rep, bool reject_nonpositive_eps) const {

  // ROOTSTAT status = ROOTSTAT::SUCCESS;
  rep.iters = 0;
  rep.adjust_cons = false;
  rep.set_atmo = false;
  rep.soft_root_conv = false;
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

  /* Copy cons vector, prevent round-off errors */
  const cons_vars cv_const = cv;

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

  if (!std::isfinite(cv.tau) || cv.dens <= 0.0) {
    rep.set_range_rho(cv.dens, 0.0);
    cv = cv_const;
    return;
  }
  const CCTK_REAL Ye_raw = cv.DYe / cv.dens;
  const CCTK_REAL Ye = std::clamp(Ye_raw, eos_3p->rgye.min, eos_3p->rgye.max);
  if (Ye != Ye_raw)
    rep.adjust_cons = true;
  cv.DYe = cv.dens * Ye;
  pv_seeds.Ye = Ye;

  /* update rho seed from cv and wlor */
  // rho consistent with cv.rho should be better guess than rho from last
  // timestep
  pv_seeds.rho = cv.dens / pv_seeds.w_lor;

  CCTK_REAL eps_last = max({pv_seeds.eps, atmo.eps_atmo});
  CCTK_REAL rho_seed = pv_seeds.rho;
  if constexpr (!std::is_same_v<EOSType, EOSX::eos_3p_idealgas>) {
    rho_seed = std::clamp(rho_seed, eos_3p->rgrho.min, eos_3p->rgrho.max);
    const auto rgeps_seed = eos_3p->range_eps_from_rho_ye(rho_seed, Ye);
    eps_last = std::clamp(pv_seeds.eps, rgeps_seed.min, rgeps_seed.max);
  }

  /* get pressure seed from updated pv_seeds.rho */
  pv_seeds.press =
      eos_3p->press_from_rho_eps_ye(rho_seed, eps_last, pv_seeds.Ye);

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
  // constexpr CCTK_REAL dw = 1. / (1. - dv);

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
    if (x[1] >= 1.0) {
      x[1] = dv;
    }
  }

  if (x[0] <= 0.0) {
    x[0] = fabs(x[0]) + Zmin;
  } else {
    if (x[0] > 1e20) {
      x[0] = x_old[0];
    }
  }

  CCTK_REAL temp = pv_seeds.temperature;
  CCTK_INT k;
  for (k = 1; k <= maxIterations; k++) {

    /* Expressions for the jacobian are adapted from the Noble C2P
    implementation in the Spritz code. As the analytical form of the equations
    is known, the Newton-Raphson step can be computed explicitly */

    const CCTK_REAL Z = x[0];
    const CCTK_REAL invZ = 1.0 / Z;
    const CCTK_REAL Vsq = x[1];

    const CCTK_REAL Sdotn = -(cv.tau + cv.dens);
    CCTK_REAL p_tmp, dPdZ, dPdvsq;
    if (!get_Press_funcZVsq(p_tmp, dPdZ, dPdvsq, Z, Vsq, eos_3p, cv, &temp)) {
      rep.set_root_conv();
      cv = cv_const;
      return;
    }

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
    if (!std::isfinite(detjac_inv)) {
      rep.set_root_conv();
      cv = cv_const;
      return;
    }

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
    CCTK_REAL step = 1.0;
    if constexpr (std::is_same_v<EOSType, EOSX::eos_3p_tabulated3d>) {
      // Keep table trials inside the EOS domain.
      bool valid = false;
      for (CCTK_INT trial = 0; trial < maxIterations; ++trial) {
        x[0] = x_old[0] + step * dx[0];
        x[1] = fmax(0.0, x_old[1] + step * dx[1]);
        CCTK_REAL temp_trial = temp;
        if (get_Press_funcZVsq(p_tmp, dPdZ, dPdvsq, x[0], x[1],
                               eos_3p, cv, &temp_trial)) {
          temp = temp_trial;
          valid = true;
          break;
        }
        step *= 0.5;
        if (step <= std::numeric_limits<CCTK_REAL>::epsilon())
          break;
      }
      if (!valid) {
        rep.set_root_conv();
        cv = cv_const;
        return;
      }
    }

    /* make sure that the new x[] is physical */
    if (x[1] < 0.0) {
      x[1] = 0.0;
    }

    else {
      if (x[1] >= 1.0) {
        x[1] = dv;
      }
    }

    if (x[0] <= 0.0) {
      x[0] = fabs(x[0]) + Zmin;
    } else {
      if (x[0] > 1e20) {
        x[0] = x_old[0];
      }
    }

    // calculate the convergence criterion
    errx = (x[0] == 0.) ? fabs(step * dx[0]) : fabs(step * dx[0] / x[0]);

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
    bool accept_soft = false;
    if (soft_root_convergence && std::isfinite(errx)) {
      const CCTK_REAL soft_tol = soft_root_width_factor * tolerance;
      accept_soft = (fabs(errx) <= soft_tol);
    }
    if (!accept_soft) {
      // set status to root not converged
      rep.set_root_conv();
      // status = ROOTSTAT::NOT_CONVERGED;
      cv = cv_const;
      return;
    }
    rep.set_soft_root_conv();
  }

  // Check for bad untrapped divergences
  if ((!isfinite(f)) || (!isfinite(df))) {
    rep.set_root_bracket();
    cv = cv_const;
    return;
  }

  /* Calculate primitives from Z and W */
  CCTK_REAL Z_Sol = x[0];
  CCTK_REAL vsq_Sol = x[1];

  if constexpr (!std::is_same_v<EOSType, EOSX::eos_3p_idealgas>) {
    // A shortened step is not proof of convergence. Check both equations
    // with the configured root tolerance and a roundoff-sized lower bound.
    CCTK_REAL press_final, dPdZ_final, dPdVsq_final;
    if (!get_Press_funcZVsq(press_final, dPdZ_final, dPdVsq_final,
                           Z_Sol, vsq_Sol, eos_3p, cv, &temp)) {
      rep.set_root_conv();
      cv = cv_const;
      return;
    }
    const CCTK_REAL bz2 = BiSi * BiSi / (Z_Sol * Z_Sol);
    const CCTK_REAL mom_resid = Ssq - vsq_Sol * pow(Bsq + Z_Sol, 2) +
                                bz2 * (Bsq + 2.0 * Z_Sol);
    const CCTK_REAL energy_resid = cv.tau + cv.dens -
        0.5 * Bsq * (1.0 + vsq_Sol) + 0.5 * bz2 - Z_Sol + press_final;
    const CCTK_REAL residual_tol = fmax(
        tolerance * (soft_root_convergence ? soft_root_width_factor : 1.0),
        32.0 * std::numeric_limits<CCTK_REAL>::epsilon());
    if (!std::isfinite(mom_resid) || !std::isfinite(energy_resid) ||
        fabs(mom_resid) > residual_tol * fmax(Ssq, pow(Bsq + Z_Sol, 2)) ||
        fabs(energy_resid) > residual_tol *
            (fabs(cv.tau) + cv.dens + Bsq + Z_Sol + fabs(press_final))) {
      rep.set_root_conv();
      cv = cv_const;
      return;
    }
  }

  /* Write prims if C2P succeeded */
  CCTK_REAL eps_raw;
  WZ2Prim(Z_Sol, vsq_Sol, Bsq, BiSi, eos_3p, pv, eps_raw, cv, gup, glo);

  // Error out if rho is negative or zero
  if (!std::isfinite(pv.rho) || pv.rho <= 0.0) {
    // set status to rho is out of range
    rep.set_range_rho(cv.dens, pv.rho);
    cv = cv_const;
    return;
  }

  // Let the usual temperature floor repair non-positive eps unless the caller
  // has an entropy-based fallback available.
  const auto rgeps = eos_3p->range_eps_from_rho_ye(pv.rho, pv.Ye);
  const bool eps_invalid = rgeps.min < 0.0
                               ? eps_raw < rgeps.min : eps_raw <= 0.0;
  if (!std::isfinite(eps_raw) || (reject_nonpositive_eps && eps_invalid)) {
    rep.set_range_eps(eps_raw);
    cv = cv_const;
    return;
  }

  if (eps_raw < rgeps.min || eps_raw > rgeps.max)
    rep.adjust_cons = true;

  // set to atmo if computed rho is below floor density
  if (pv.rho < atmo.rho_cut) {
    rep.set_atmo_set();
    atmo.set(pv, cv, glo);
    return;
  }

  c2p::prims_floors_and_ceilings(eos_3p, pv, cv, alp, beta, glo, rep);
  if (rep.failed()) {
    cv = cv_const;
    return;
  }

  // Recompute cons if prims have been adjusted
  if (rep.adjust_cons) {
    cv.from_prim(pv, glo);
    cv.dBvec = cv_const.dBvec;
  } else {
    cv = cv_const;
    // Conserved entropy must be consistent with new prims
    cv.DEnt = cv.dens * pv.entropy;
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
