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

#include "eos.hxx"
#include "eos_idealgas.hxx"

#include "c2p_utils.hxx"

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

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
      get_Ssq_Exact(const vec<CCTK_REAL, 3> &mom,
                    const smat<CCTK_REAL, 3> &gup) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
      get_Bsq_Exact(const vec<CCTK_REAL, 3> &B_up,
                    const smat<CCTK_REAL, 3> &glo) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
      get_BiSi_Exact(const vec<CCTK_REAL, 3> &Bvec, const vec<CCTK_REAL, 3> &mom) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 2>
      get_WLorentz_bsq_Seeds(const vec<CCTK_REAL, 3> &B_up, const vec<CCTK_REAL, 3> &v_up,
                             const smat<CCTK_REAL, 3> &glo) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  prims_floors_and_ceilings(const EOSType &eos_th, prim_vars &pv, const cons_vars &cv,
        const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const;

  public:

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  bh_interior_fail(const EOSType &eos_th, prim_vars &pv, cons_vars &cv,
        const smat<CCTK_REAL, 3> &glo) const;
};

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p::prims_floors_and_ceilings(const EOSType &eos_th, prim_vars &pv, const cons_vars &cv,
                               const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const {

  //const CCTK_REAL spatial_detg = calc_det(glo);
  //const CCTK_REAL sqrt_detg = sqrt(spatial_detg);

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
    pv.eps = eos_th.eps_from_valid_rho_press_ye(pv.rho, pv.press, pv.Ye);
    pv.entropy =
        eos_th.kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
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

  if (pv.rho > eos_th.rgrho.max) {

    // remove mass, changes conserved density D
    pv.rho = eos_th.rgrho.max;
    pv.eps = eos_th.eps_from_valid_rho_press_ye(eos_th.rgrho.max, pv.press, pv.Ye);
    pv.entropy =
        eos_th.kappa_from_valid_rho_eps_ye(eos_th.rgrho.max, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  }

  // ----------
  // Floor and ceiling for eps
  // Keeps rho the same and changes press
  // ----------

  // check the validity of the computed eps
  auto rgeps = eos_th.range_eps_from_valid_rho_ye(pv.rho, pv.Ye);
  if (pv.eps > rgeps.max) {
    // printf("(pv.eps > rgeps.max) is true, adjusting cons.. \n");
    // if (pv.rho >= rho_strict) {
    //  rep.set_range_eps(pv.eps); // sets adjust_cons to false by default
    //  rep.adjust_cons = true;
    //  set_to_nan(pv, cv);
    //  return;
    //}
    pv.eps = rgeps.max;
    pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.entropy = eos_th.kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  } else if (pv.eps < rgeps.min) {
    /*
    printf(
        "(pv.eps < rgeps.min) is true! pv.eps, rgeps.min: %26.16e, %26.16e
    \n",
        pv.eps, rgeps.min);
    printf(" Not adjusting cons.. \n");
    */
    // rep.set_range_eps(rgeps.min); // sets adjust_cons to true
    pv.eps = rgeps.min;
    pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.entropy = eos_th.kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  }

  // TODO: check validity for Ye

}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p::bh_interior_fail(const EOSType &eos_th, prim_vars &pv, cons_vars &cv,
                      const smat<CCTK_REAL, 3> &glo) const {

  // Treatment for BH interiors after C2P failures
  // NOTE: By default, alp_thresh=0 so the if condition below is never
  // triggered. One must be very careful when using this functionality and
  // must correctly set alp_thresh, rho_BH, eps_BH and vwlim_BH in the
  // parfile
  pv.rho = rho_BH; // typically set to 0.01% to 1% of rho_max of initial
                   // NS or disk
  pv.eps = eps_BH;
  pv.Ye = 0.5;
  pv.press =
      eos_th.press_from_valid_rho_eps_ye(rho_BH, eps_BH, 0.5);
  pv.entropy =
      eos_th.kappa_from_valid_rho_eps_ye(rho_BH, eps_BH, 0.5);
  // check on velocities
  CCTK_REAL wlim_BH = sqrt(1.0 + vwlim_BH * vwlim_BH);
  CCTK_REAL vlim_BH = vwlim_BH / wlim_BH;
  CCTK_REAL sol_v = sqrt((pv.w_lor * pv.w_lor - 1.0)) / pv.w_lor;
  if (sol_v > vlim_BH) {
    pv.vel *= vlim_BH / sol_v;
    pv.w_lor = wlim_BH;
  }
  cv.from_prim(pv, glo);

}

} // namespace Con2PrimFactory

#endif
