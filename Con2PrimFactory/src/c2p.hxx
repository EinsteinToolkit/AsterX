/*! \file c2p.hxx
\brief Defines a c2p
\author Jay Kalinani

c2p is effectively an interface to be used by different c2p implementations.

*/

#ifndef C2P_HXX
#define C2P_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <math.h>

#include "prims.hxx"
#include "cons.hxx"
#include "atmo.hxx"
#include "c2p_report.hxx"
#include "c2p_utils.hxx"
#include "setup_eos.hxx"

namespace Con2PrimFactory {

using namespace AsterUtils;

constexpr CCTK_INT X = 0;
constexpr CCTK_INT Y = 1;
constexpr CCTK_INT Z = 2;

/* Abstract class c2p */
class c2p {
protected:
  /* The constructor must initialize the following variables */

  atmosphere atmo;
  CCTK_INT maxIterations;
  CCTK_REAL tolerance;
  CCTK_REAL alp_thresh;
  CCTK_REAL cons_error;
  CCTK_REAL vw_lim;
  CCTK_REAL w_lim;
  CCTK_REAL v_lim;
  CCTK_REAL Bsq_lim;
  CCTK_REAL rho_BH;
  CCTK_REAL eps_BH;
  CCTK_REAL vwlim_BH;
  bool ye_lenient;
  bool use_zprim;
  bool use_temp;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Ssq_Exact(const vec<CCTK_REAL, 3> &mom,
                const smat<CCTK_REAL, 3> &gup) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Bsq_Exact(const vec<CCTK_REAL, 3> &B_up,
                const smat<CCTK_REAL, 3> &glo) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_BiSi_Exact(const vec<CCTK_REAL, 3> &Bvec,
                 const vec<CCTK_REAL, 3> &mom) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 2>
  get_WLorentz_bsq_Seeds(const vec<CCTK_REAL, 3> &B_up,
                         const vec<CCTK_REAL, 3> &v_up,
                         const smat<CCTK_REAL, 3> &glo) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  prims_floors_and_ceilings(const EOSType *eos_3p, prim_vars &pv,
                            const cons_vars &cv, const smat<CCTK_REAL, 3> &glo,
                            c2p_report &rep) const;

public:
  template <typename EOSType, bool limiting>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  bh_interior(const EOSType *eos_3p, prim_vars &pv, cons_vars &cv,
              const smat<CCTK_REAL, 3> &glo) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  cons_floors_and_ceilings(const EOSType *eos_3p, cons_vars &cv, 
                           const smat<CCTK_REAL, 3> &glo) const;
};

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p::prims_floors_and_ceilings(const EOSType *eos_3p, prim_vars &pv,
                               const cons_vars &cv,
                               const smat<CCTK_REAL, 3> &glo,
                               c2p_report &rep) const {

  // const CCTK_REAL spatial_detg = calc_det(glo);
  // const CCTK_REAL sqrt_detg = sqrt(spatial_detg);

  // Lower velocity
  const vec<CCTK_REAL, 3> v_low = calc_contraction(glo, pv.vel);

  // ----------
  // Floor and ceiling for rho and velocity
  // Keeps pressure the same and changes eps
  // ----------

  // check if computed velocities are within the specified limit
  CCTK_REAL vsq_Sol = calc_contraction(v_low, pv.vel);
  CCTK_REAL sol_v = sqrt(vsq_Sol);
  if (sol_v > v_lim) {
    /*
    printf("(sol_v > v_lim) is true! \n");
    printf("sol_v, v_lim: %26.16e, %26.16e \n", sol_v, v_lim);
    */
    // add mass, keeps conserved density D
    pv.rho = cv.dens / w_lim;
    if (use_temp) {
      pv.eps =
          eos_3p->eps_from_valid_rho_temp_ye(pv.rho, pv.temperature, pv.Ye);
      pv.press =
          eos_3p->press_from_valid_rho_temp_ye(pv.rho, pv.temperature, pv.Ye);
    } else {
      pv.eps = eos_3p->eps_from_valid_rho_press_ye(pv.rho, pv.press, pv.Ye);
    }
    pv.entropy = eos_3p->kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    // if (pv.rho >= rho_strict) {
    //  rep.set_speed_limit({ sol_v, sol_v, sol_v });
    //  set_to_nan(pv, cv);
    //  return;
    //}
    pv.vel *= v_lim / sol_v;
    pv.w_lor = w_lim;
    // pv.eps = std::min(std::max(eos_th.rgeps.min, pv.eps),
    // eos_th.rgeps.max);
    // pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  }

  if (pv.rho > eos_3p->rgrho.max) {

    // remove mass, changes conserved density D
    pv.rho = eos_3p->rgrho.max;
    if (use_temp) {
      pv.eps =
          eos_3p->eps_from_valid_rho_temp_ye(pv.rho, pv.temperature, pv.Ye);
      pv.press =
          eos_3p->press_from_valid_rho_temp_ye(pv.rho, pv.temperature, pv.Ye);
    } else {
      pv.eps = eos_3p->eps_from_valid_rho_press_ye(pv.rho, pv.press, pv.Ye);
    }
    pv.entropy =
        eos_3p->kappa_from_valid_rho_eps_ye(eos_3p->rgrho.max, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  }

  // ----------
  // Floor and ceiling for eps
  // Keeps rho the same and changes press
  // ----------

  // check the validity of the computed eps
  auto rgeps = eos_3p->range_eps_from_valid_rho_ye(pv.rho, pv.Ye);
  
  if (pv.eps > rgeps.max) {
    // printf("(pv.eps > rgeps.max) is true, adjusting cons.. \n");
    // if (pv.rho >= rho_strict) {
    //  rep.set_range_eps(pv.eps); // sets adjust_cons to false by default
    //  rep.adjust_cons = true;
    //  set_to_nan(pv, cv);
    //  return;
    //}
    pv.eps = rgeps.max;
    pv.temperature = eos_3p->temp_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.press = eos_3p->press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.entropy = eos_3p->kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  } else if (pv.eps < atmo.eps_atmo) {
    /*
    printf(
        "(pv.eps < rgeps.min) is true! pv.eps, rgeps.min: %26.16e, %26.16e
    \n",
        pv.eps, rgeps.min);
    printf(" Not adjusting cons.. \n");
    */
    // rep.set_range_eps(rgeps.min); // sets adjust_cons to true
    pv.eps = atmo.eps_atmo;
    pv.temperature = eos_3p->temp_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.press = eos_3p->press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.entropy = eos_3p->kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  }

  // TODO: check validity for Ye
}

template <typename EOSType, bool limiting>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p::bh_interior(const EOSType *eos_3p, prim_vars &pv, cons_vars &cv,
                 const smat<CCTK_REAL, 3> &glo) const {

  // Treatment for BH interiors after C2P failures
  // NOTE: By default, alp_thresh=0 so the if condition below is never
  // triggered. One must be very careful when using this functionality and
  // must correctly set alp_thresh, rho_BH, eps_BH and vwlim_BH in the
  // parfile

  const CCTK_REAL wlim_BH = sqrt(1.0 + vwlim_BH * vwlim_BH);
  const CCTK_REAL vlim_BH = vwlim_BH / wlim_BH;

  bool recomp_flag = false;

  if constexpr (limiting) {

    if (pv.rho > rho_BH) {
      pv.rho = rho_BH; // typically set to 0.01% to 1% of rho_max of initial
                       // NS or disk
      recomp_flag = true;
    };

    if (pv.eps > eps_BH) {
      pv.eps = eps_BH;
      recomp_flag = true;
    };

    const CCTK_REAL sol_v = sqrt((pv.w_lor * pv.w_lor - 1.0)) / pv.w_lor;
    if (sol_v > vlim_BH) {
      pv.vel *= vlim_BH / sol_v;
      pv.w_lor = wlim_BH;
      recomp_flag = true;
    };

    if (recomp_flag) {

      pv.Ye = atmo.ye_atmo;
  
      pv.temperature = eos_3p->temp_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
      pv.press = eos_3p->press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye); 
      pv.entropy = eos_3p->kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
  
      cv.from_prim(pv, glo);
    };
  
  } else {
 
    pv.rho = rho_BH; // typically set to 0.01% to 1% of rho_max of initial
                     // NS or disk
    pv.eps = eps_BH;
    pv.Ye = atmo.ye_atmo;
  
    pv.temperature = eos_3p->temp_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.press = eos_3p->press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye); 
    pv.entropy = eos_3p->kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    // Set velocity such that new conserved momentum has same 
    // direction as before

    // Inverse metric
    const CCTK_REAL spatial_detg = calc_det(glo);
    const smat<CCTK_REAL, 3> gup = calc_inv(glo, spatial_detg);

    // Compute Z = rho * h * W * W
    const CCTK_REAL Z_loc = ( pv.rho * ( 1.0 + pv.eps ) + pv.press ) * wlim_BH * wlim_BH;

    // Get Bsq
    const vec<CCTK_REAL, 3> B_low = calc_contraction(glo, pv.Bvec);
    const CCTK_REAL Bsq = calc_contraction(B_low, pv.Bvec);

    // Norm of conserved momentum, undensitize here
    vec<CCTK_REAL, 3> mom_low = cv.mom / sqrt(spatial_detg);
    vec<CCTK_REAL, 3> mom_up  = calc_contraction(gup, mom_low);
    const CCTK_REAL Ssq_old = calc_contraction(mom_low, mom_up);
    const CCTK_REAL S_old = sqrt(Ssq_old) + 1e-50;

    // Get BiSi = S_iB^i
    const CCTK_REAL BiSi_old = calc_contraction(mom_low, pv.Bvec);

    // Normalize S_iB^i by S = sqrt(S_iS^i)
    const CCTK_REAL BiEsi = BiSi_old / S_old;

    // Compute magnitude of new conserved momentum
    const CCTK_REAL Ssq_new = ( (Z_loc + Bsq)*(Z_loc + Bsq)*vlim_BH*vlim_BH ) / 
                              ( 1.0 + BiEsi * BiEsi * ( 2.0 * Z_loc + Bsq ) / ( Z_loc * Z_loc ) );
    const CCTK_REAL S_new = sqrt(Ssq_new);

    // Rescale momenta
    mom_low *= S_new / S_old;
    mom_up  *= S_new / S_old;

    // Finally, compute velocity 
    // This is (24) from https://arxiv.org/pdf/1712.07538
    pv.vel(X) = mom_up(X) /
                (Z_loc + Bsq);
    pv.vel(X) += BiEsi * S_new * pv.Bvec(X) / (Z_loc * (Z_loc + Bsq));

    pv.vel(Y) = mom_up(Y) /
                (Z_loc + Bsq);
    pv.vel(Y) += BiEsi * S_new * pv.Bvec(Y) / (Z_loc * (Z_loc + Bsq));

    pv.vel(Z) = mom_up(Z) /
                (Z_loc + Bsq);
    pv.vel(Z) += BiEsi * S_new * pv.Bvec(Z) / (Z_loc * (Z_loc + Bsq));

    pv.w_lor = wlim_BH;

    cv.from_prim(pv, glo);
 
  };
};

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p::cons_floors_and_ceilings(const EOSType *eos_3p, cons_vars &cv, 
                              const smat<CCTK_REAL, 3> &glo) const {

  // Limit conservative variables
  // Note that conservatives are densitized

  const CCTK_REAL spatial_detg = calc_det(glo);
  const CCTK_REAL sqrt_detg = sqrt(spatial_detg);

  const smat<CCTK_REAL, 3> gup = calc_inv(glo, spatial_detg);

  // Lower limit on tau/conserved internal energy
  // Based on Appendix A of https://arxiv.org/pdf/1112.0568

  // Compute Bsq
  vec<CCTK_REAL, 3> B_low  = calc_contraction(glo, cv.dBvec);
  const CCTK_REAL BsqL = calc_contraction(B_low, cv.dBvec);
  const CCTK_REAL tau_lim = 0.5*BsqL/sqrt_detg;

  if (cv.tau <= tau_lim) {
    cv.tau = tau_lim;
  }

  // Dominant energy condition 
  // (A5) from https://arxiv.org/pdf/1505.01607

  vec<CCTK_REAL, 3> mom_up  = calc_contraction(gup, cv.mom);
  const CCTK_REAL mom2L = calc_contraction(cv.mom, mom_up);

  const CCTK_REAL slim  = cv.dens + cv.tau;
  const CCTK_REAL slim2 = slim*slim;

  if (mom2L > slim2) {
   // (A51) from https://arxiv.org/pdf/1112.0568
   cv.mom = cv.mom * sqrt(slim2/mom2L);
  }

};

} // namespace Con2PrimFactory

#endif
