#ifndef ASTERX_RECONSTRUCT_HXX
#define ASTERX_RECONSTRUCT_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "utils.hxx"

namespace AsterX {

using namespace std;
using namespace Loop;

// Struct used to pass parameters to the reconstruction routine
typedef struct {
  // PPM parameters
  bool ppm_shock_detection, ppm_zone_flattening;
  CCTK_REAL poly_k, poly_gamma;
  CCTK_REAL ppm_eta1, ppm_eta2;
  CCTK_REAL ppm_eps;
  CCTK_REAL ppm_omega1, ppm_omega2;
  // WENOZ parameters
  CCTK_REAL weno_eps;
} reconstruct_params_t;

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T minmod(const T &x,
                                                                   const T &y) {
  if (signbit(x) != signbit(y))
    return T(0);
  if (fabs(x) < fabs(y))
    return x;
  else
    return y;
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
monocentral(const T &x, const T &y) {
  if (sgn(x) != sgn(y))
    return 0;
  else
    return sgn(x) * min(2 * fabs(x), min(2 * fabs(y), fabs(x + y) / 2));
}

/* PPM reconstruction scheme. See Colella & Woodward (1984) (e.g. at
 * https://crd.lbl.gov/assets/pubs_presos/AMCS/ANAG/A141984.pdf) */
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<CCTK_REAL, 2>
ppm(const GF3D2<const CCTK_REAL> &gf_var,
    const array<const vect<int, dim>, 5> &cells, const CCTK_INT &dir,
    const bool &gf_is_rho, const GF3D2<const CCTK_REAL> &gf_press,
    const GF3D2<const CCTK_REAL> &gf_vel_dir,
    const reconstruct_params_t &reconstruct_params) {
  // Unpack all cells in the stencil
  const auto &Imm = cells.at(0);
  const auto &Im = cells.at(1);
  const auto &I = cells.at(2);
  const auto &Ip = cells.at(3);
  const auto &Ipp = cells.at(4);

  // Unpack all PPM parameters
  const bool &ppm_shock_detection = reconstruct_params.ppm_shock_detection;
  const bool &ppm_zone_flattening = reconstruct_params.ppm_zone_flattening;
  const CCTK_REAL &poly_k = reconstruct_params.poly_k;
  const CCTK_REAL &poly_gamma = reconstruct_params.poly_gamma;
  const CCTK_REAL &ppm_eta1 = reconstruct_params.ppm_eta1;
  const CCTK_REAL &ppm_eta2 = reconstruct_params.ppm_eta2;
  const CCTK_REAL &ppm_eps = reconstruct_params.ppm_eps;
  const CCTK_REAL &ppm_omega1 = reconstruct_params.ppm_omega1;
  const CCTK_REAL &ppm_omega2 = reconstruct_params.ppm_omega2;

  // Grid function at neighboring cells
  const CCTK_REAL &gf_Imm = gf_var(Imm);
  const CCTK_REAL &gf_Im = gf_var(Im);
  const CCTK_REAL &gf_I = gf_var(I);
  const CCTK_REAL &gf_Ip = gf_var(Ip);
  const CCTK_REAL &gf_Ipp = gf_var(Ipp);

  // Helpers
  const CCTK_REAL diff_Im = gf_I - gf_Imm;
  const CCTK_REAL diff_I = gf_Ip - gf_Im;
  const CCTK_REAL diff_Ip = gf_Ipp - gf_I;

  const CCTK_REAL delta_Im = 0.5 * diff_Im;
  const CCTK_REAL delta_I = 0.5 * diff_I;
  const CCTK_REAL delta_Ip = 0.5 * diff_Ip;

  const CCTK_REAL _2fabs_Im_Imm = 2 * fabs(gf_Im - gf_Imm);
  const CCTK_REAL _2fabs_I_Im = 2 * fabs(gf_I - gf_Im);
  const CCTK_REAL _2fabs_Ip_I = 2 * fabs(gf_Ip - gf_I);
  const CCTK_REAL _2fabs_Ipp_Ip = 2 * fabs(gf_Ipp - gf_Ip);

  const bool same_sgn_Im = ((gf_I - gf_Im) * (gf_Im - gf_Imm) > 0);
  const bool same_sgn_I = ((gf_Ip - gf_I) * (gf_I - gf_Im) > 0);
  const bool same_sgn_Ip = ((gf_Ipp - gf_Ip) * (gf_Ip - gf_I) > 0);

  const CCTK_REAL deltamod_Im =
      same_sgn_Im
          ? sgn(delta_Im) * min(fabs(delta_Im), min(_2fabs_Im_Imm, _2fabs_I_Im))
          : 0;
  const CCTK_REAL deltamod_I =
      same_sgn_I
          ? sgn(delta_I) * min(fabs(delta_I), min(_2fabs_I_Im, _2fabs_Ip_I))
          : 0;
  const CCTK_REAL deltamod_Ip =
      same_sgn_Ip
          ? sgn(delta_Ip) * min(fabs(delta_Ip), min(_2fabs_Ip_I, _2fabs_Ipp_Ip))
          : 0;

  /* Initial reconstructed states at the interfaces between cells Im/I ans I/Ip
   * NOTE: not const because they may change later */
  const CCTK_REAL gf_Im_I =
      0.5 * (gf_Im + gf_I) + (1. / 6.) * (deltamod_Im - deltamod_I);
  const CCTK_REAL gf_I_Ip =
      0.5 * (gf_I + gf_Ip) + (1. / 6.) * (deltamod_I - deltamod_Ip);

  CCTK_REAL rc_low = gf_Im_I;
  CCTK_REAL rc_up = gf_I_Ip;

  /* Pressure may or may not be needed, or may not be needed at all points
   * (depending on whether shock detection and zone flattening are activated or
   * not) */
  // const CCTK_REAL &press_Immm = gf_press(Immm);  // Only used in the original
  // zone flattening scheme
  const CCTK_REAL &press_Imm = gf_press(Imm);
  const CCTK_REAL &press_Im = gf_press(Im);
  // const CCTK_REAL &press_I    = gf_press(I);     // Only used in the original
  // zone flattening scheme
  const CCTK_REAL &press_Ip = gf_press(Ip);
  const CCTK_REAL &press_Ipp = gf_press(Ipp);
  // const CCTK_REAL &press_Ippp = gf_press(Ippp);  // Only used in the original
  // zone flattening scheme

  const CCTK_REAL diff_press_I = press_Ip - press_Im;
  const CCTK_REAL min_press_I = min(press_Im, press_Ip);

  /* Shock detection (eqs. 1.15, 1.16, 1.17 with uniform grid spacing).
   * This is only applied to rho and only if the shock is marked as a contact
   * discontinuity (see eq. 3.2) */
  // FIXME: contact discontinuity check only valid for polytropic/ideal-fluid
  // EOS
  if (ppm_shock_detection and gf_is_rho) {
    const CCTK_REAL k_gamma = poly_k * poly_gamma;
    const bool contact_disc = (k_gamma * fabs(diff_I) * min_press_I >=
                               fabs(diff_press_I) * min(gf_Ip, gf_Im));
    if (contact_disc) {
      // This assumes gf_var is rho (gf_is_rho is true)
      const CCTK_REAL d2rho_Im = gf_I - 2 * gf_Im + gf_Imm;
      const CCTK_REAL d2rho_Ip = gf_Ipp - 2 * gf_Ip + gf_I;
      const CCTK_REAL d2_prod = d2rho_Im * d2rho_Ip;
      const CCTK_REAL eta_tilde_I =
          (d2_prod < 0) ? (-1. / 6.) * d2_prod / diff_I : 0;
      const CCTK_REAL eta_I =
          max(0., min(ppm_eta1 * (eta_tilde_I - ppm_eta2), 1.));

      rc_low = (1 - eta_I) * rc_low + eta_I * (gf_Im + 0.5 * deltamod_Im);
      rc_up = (1 - eta_I) * rc_up + eta_I * (gf_Ip - 0.5 * deltamod_Ip);
    }
  }

  /* Zone flattening to reduce post-shock oscillations (see eq. 4.1 + appendix).
   * The flattening parameter f_I is set to ftilde_I instead of
   * max(ftilde_I, ftilde_{I+s_I}) (where s_I can be +1 or -1), thereby avoiding
   * using four ghost cells at interprocess/domain boundaries. This should not
   * be a major issue and is done in many GRMHD codes (e.g. WhiskyMHD, GRHydro,
   * Spritz, IllinoisGRMHD). */
  if (ppm_zone_flattening) {
    const CCTK_REAL w_I = (diff_press_I > ppm_eps * min_press_I and
                           gf_vel_dir(Im) > gf_vel_dir(Ip))
                              ? 1
                              : 0;
    const CCTK_REAL ftilde_I =
        1 - max(0., w_I * ppm_omega2 *
                        (diff_press_I / (press_Ipp - press_Imm) - ppm_omega1));

    const CCTK_REAL one_minus_ftilde_I = 1 - ftilde_I;
    const CCTK_REAL ftildeI_gfI = ftilde_I * gf_I;

    rc_low = ftildeI_gfI + one_minus_ftilde_I * rc_low;
    rc_up = ftildeI_gfI + one_minus_ftilde_I * rc_up;

    // This would require one more ghost cell and it's not worth it
    /*if (diff_press_I < 0) {
      const CCTK_REAL diff_press_Ip = press_Ipp - press_I;
      const CCTK_REAL w_Ip          = (diff_press_Ip > ppm_eps*min(press_I,
    press_Ipp) and gf_vel_dir(I) > gf_vel_dir(Ipp)) ? 1 : 0; const CCTK_REAL
    ftilde_Ip     = 1 - max(0., w_Ip*ppm_omega2*(diff_press_Ip/(press_Ippp -
    press_Im) - ppm_omega1)); const CCTK_REAL f_I           = max(ftilde_I,
    ftilde_Ip); const CCTK_REAL one_minus_fI  = 1 - f_I; const CCTK_REAL fI_gfI
    = f_I*gf_I; rc_low = fI_gfI + one_minus_fI*rc_low; rc_up  = fI_gfI +
    one_minus_fI*rc_up;
    }

    else {
      const CCTK_REAL diff_press_Im = press_I - press_Imm;
      const CCTK_REAL w_Im          = (diff_press_Im > ppm_eps*min(press_Imm,
    press_I) and gf_vel_dir(Imm) > gf_vel_dir(I)) ? 1 : 0; const CCTK_REAL
    ftilde_Im     = 1 - max(0., w_Im*ppm_omega2*(diff_press_Im/(press_Ip -
    press_Immm) - ppm_omega1)); const CCTK_REAL f_I           = max(ftilde_I,
    ftilde_Im); const CCTK_REAL one_minus_fI  = 1 - f_I; const CCTK_REAL fI_gfI
    = f_I*gf_I; rc_low = fI_gfI + one_minus_fI*rc_low; rc_up  = fI_gfI +
    one_minus_fI*rc_up;
    }*/
  }

  // Monotonization (see eq. 1.10)
  if ((rc_up - gf_I) * (gf_I - rc_low) <= 0) {
    rc_low = gf_I;
    rc_up = gf_I;
  }

  else {
    const CCTK_REAL diff_rc = rc_up - rc_low;
    const CCTK_REAL diff_rc_sq = diff_rc * diff_rc;
    const CCTK_REAL gf6_I = 6 * diff_rc * (gf_I - 0.5 * (rc_low + rc_up));

    if (gf6_I > diff_rc_sq) {
      rc_low = 3 * gf_I - 2 * rc_up;
    }

    else if (gf6_I < -diff_rc_sq) {
      rc_up = 3 * gf_I - 2 * rc_low;
    }
  }

  // Return the lower and upper reconstructed states in cell I
  const array<CCTK_REAL, 2> rc = {rc_low, rc_up};
  return rc;
}

/* WENO-Z: Performs the reconstruction of a given variable using 5th order
 * WENO-Z method based on Borges et al., ''An improved weighted essentially
 * non-oscillatory scheme for hyperbolic conservation laws'', 2008) */
// Also, see the Spritz code for details

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<CCTK_REAL, 2>
wenoz(const GF3D2<const CCTK_REAL> &gf_var,
      const array<const vect<int, dim>, 5> &cells,
      const reconstruct_params_t &reconstruct_params) {

  // Unpack all WENO-Z parameters
  const CCTK_REAL weno_eps = reconstruct_params.weno_eps;

  // Unpack all cells in the stencil
  const auto &Imm = cells.at(0);
  const auto &Im = cells.at(1);
  const auto &I = cells.at(2);
  const auto &Ip = cells.at(3);
  const auto &Ipp = cells.at(4);

  // Grid function at neighboring cells
  const CCTK_REAL &gf_Imm = gf_var(Imm);
  const CCTK_REAL &gf_Im = gf_var(Im);
  const CCTK_REAL &gf_I = gf_var(I);
  const CCTK_REAL &gf_Ip = gf_var(Ip);
  const CCTK_REAL &gf_Ipp = gf_var(Ipp);

  // Computing the smoothness indicators (Borges et al. 2008)
  const vec<CCTK_REAL, 3> betaZ{
      (13.0 / 12.0) * pow2(gf_Imm - 2.0 * gf_Im + gf_I) +
          (1.0 / 4.0) * pow2(gf_Imm - 4.0 * gf_Im + 3.0 * gf_I),

      (13.0 / 12.0) * pow2(gf_Im - 2.0 * gf_I + gf_Ip) +
          (1.0 / 4.0) * pow2(gf_Im - gf_Ip),

      (13.0 / 12.0) * pow2(gf_I - 2.0 * gf_Ip + gf_Ipp) +
          (1.0 / 4.0) * pow2(3.0 * gf_I - 4.0 * gf_Ip + gf_Ipp)};

  // Defining tau5 based on eq. 25 of (Borges et al. 2008)
  const CCTK_REAL tau5 = abs(betaZ(0) - betaZ(2));

  // Unnormalized weights
  // Optimal weights are chosen according to Del Zanna et al. 2007

  vec<CCTK_REAL, 3> aux_alphaZ;
  aux_alphaZ(0) = 1.0 + tau5 / (betaZ(0) + weno_eps);
  aux_alphaZ(1) = 1.0 + tau5 / (betaZ(1) + weno_eps);
  aux_alphaZ(2) = 1.0 + tau5 / (betaZ(2) + weno_eps);
  const vec<CCTK_REAL, 3> wt = {5.0 / 16.0, 10.0 / 16.0, 1.0 / 16.0};
  vec<vec<CCTK_REAL, 2>, 3> alphaZ;

  // for minus side
  alphaZ(0)(0) = wt(0) * aux_alphaZ(0);
  alphaZ(1)(0) = wt(1) * aux_alphaZ(1);
  alphaZ(2)(0) = wt(2) * aux_alphaZ(2);

  // for plus side
  alphaZ(0)(1) = wt(2) * aux_alphaZ(0);
  alphaZ(1)(1) = wt(1) * aux_alphaZ(1);
  alphaZ(2)(1) = wt(0) * aux_alphaZ(2);

  // Normalized weights for the reconstruction (Del Zanna et al. 2007)
  const vec<CCTK_REAL, 2> omega_denom([&](int j) ARITH_INLINE {
    return alphaZ(0)(j) + alphaZ(1)(j) + alphaZ(2)(j);
  });

  const vec<vec<CCTK_REAL, 2>, 3> omegaZ([&](int j) ARITH_INLINE {
    return vec<CCTK_REAL, 2>(
        [&](int f) ARITH_INLINE { return alphaZ(j)(f) / omega_denom(f); });
  });

  // Reconstruct cell-centered variable to left (minus) and right (plus) cell
  // interfaces

  const CCTK_REAL var_m =
      (omegaZ(2)(0) / 8.0) * (3.0 * gf_Ipp - 10.0 * gf_Ip + 15.0 * gf_I) +
      (omegaZ(1)(0) / 8.0) * (-1.0 * gf_Ip + 6.0 * gf_I + 3.0 * gf_Im) +
      (omegaZ(0)(0) / 8.0) * (3.0 * gf_I + 6.0 * gf_Im - 1.0 * gf_Imm);

  const CCTK_REAL var_p =
      (omegaZ(0)(1) / 8.0) * (3.0 * gf_Imm - 10.0 * gf_Im + 15.0 * gf_I) +
      (omegaZ(1)(1) / 8.0) * (-1.0 * gf_Im + 6.0 * gf_I + 3.0 * gf_Ip) +
      (omegaZ(2)(1) / 8.0) * (3.0 * gf_I + 6.0 * gf_Ip - 1.0 * gf_Ipp);

  const array<CCTK_REAL, 2> var_rc = {var_m, var_p};
  return var_rc;
}

// enum class for different reconstruction routines

enum class reconstruction_t { Godunov, minmod, monocentral, ppm, wenoz };

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE array<CCTK_REAL, 2>
reconstruct(const GF3D2<const CCTK_REAL> &gf_var, const PointDesc &p,
            const reconstruction_t &reconstruction, const int &dir,
            const bool &gf_is_rho, const GF3D2<const CCTK_REAL> &gf_press,
            const GF3D2<const CCTK_REAL> &gf_vel_dir,
            const reconstruct_params_t &reconstruct_params) {
  // Neighbouring "plus" and "minus" cell indices
  const auto Immm = p.I - 3 * p.DI[dir];
  const auto Imm = p.I - 2 * p.DI[dir];
  const auto Im = p.I - p.DI[dir];
  const auto Ip = p.I;
  const auto Ipp = p.I + p.DI[dir];
  const auto Ippp = p.I + 2 * p.DI[dir];

  switch (reconstruction) {

  case reconstruction_t::Godunov: {
    CCTK_REAL var_m = gf_var(Im);
    CCTK_REAL var_p = gf_var(Ip);
    return array<CCTK_REAL, 2>{var_m, var_p};
  }

  case reconstruction_t::minmod: {
    // reconstructs values of Im and Ip at the common face between these
    // two cells
    CCTK_REAL var_slope_p = gf_var(Ipp) - gf_var(Ip);
    CCTK_REAL var_slope_c = gf_var(Ip) - gf_var(Im);
    CCTK_REAL var_slope_m = gf_var(Im) - gf_var(Imm);
    // reconstructed Im on its "plus/right" side
    CCTK_REAL var_m = gf_var(Im) + minmod(var_slope_c, var_slope_m) / 2;
    // reconstructed Ip on its "minus/left" side
    CCTK_REAL var_p = gf_var(Ip) - minmod(var_slope_p, var_slope_c) / 2;
    return array<CCTK_REAL, 2>{var_m, var_p};
  }

  case reconstruction_t::monocentral: {
    // reconstructs values of Im and Ip at the common face between these
    // two cells
    // reconstructed Im on its "plus/right" side
    CCTK_REAL var_slope_p = gf_var(Ip) - gf_var(Im);
    CCTK_REAL var_slope_m = gf_var(Im) - gf_var(Imm);
    CCTK_REAL var_m = gf_var(Im) + monocentral(var_slope_p, var_slope_m) / 2;
    // reconstructed Ip on its "minus/left" side
    var_slope_p = gf_var(Ipp) - gf_var(Ip);
    var_slope_m = gf_var(Ip) - gf_var(Im);
    CCTK_REAL var_p = gf_var(Ip) - monocentral(var_slope_p, var_slope_m) / 2;
    return array<CCTK_REAL, 2>{var_m, var_p};
  }

  case reconstruction_t::ppm: {
    const array<const vect<int, dim>, 5> cells_Im = {Immm, Imm, Im, Ip, Ipp};
    const array<const vect<int, dim>, 5> cells_Ip = {Imm, Im, Ip, Ipp, Ippp};

    const array<CCTK_REAL, 2> rc_Im =
        ppm(gf_var, cells_Im, dir, gf_is_rho, gf_press, gf_vel_dir,
            reconstruct_params);
    const array<CCTK_REAL, 2> rc_Ip =
        ppm(gf_var, cells_Ip, dir, gf_is_rho, gf_press, gf_vel_dir,
            reconstruct_params);

    return array<CCTK_REAL, 2>{rc_Im.at(1), rc_Ip.at(0)};
  }

  case reconstruction_t::wenoz: {
    const array<const vect<int, dim>, 5> cells_Im = {Immm, Imm, Im, Ip, Ipp};
    const array<const vect<int, dim>, 5> cells_Ip = {Imm, Im, Ip, Ipp, Ippp};

    const array<CCTK_REAL, 2> rc_Im =
        wenoz(gf_var, cells_Im, reconstruct_params);
    const array<CCTK_REAL, 2> rc_Ip =
        wenoz(gf_var, cells_Ip, reconstruct_params);

    return array<CCTK_REAL, 2>{rc_Im.at(1), rc_Ip.at(0)};
  }

  default:
    assert(0);
  }
}

} // namespace AsterX

#endif // ASTERX_RECONSTRUCT_HXX
