#ifndef RECONX_PPM_HXX
#define RECONX_PPM_HXX

#include <cctk.h>

#include "reconx_utils.hxx"

#include <array>
#include <cmath>

namespace ReconX {

using std::array;

/**
 * @brief PPM reconstruction scheme. See Colella & Woodward (1984) (e.g. at
 * https://crd.lbl.gov/assets/pubs_presos/AMCS/ANAG/A141984.pdf)
 *
 */

template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<T, 2>
ppm(T gf_Imm, T gf_Im, T gf_I, T gf_Ip, T gf_Ipp, T press_Imm, T press_Im,
    T press_Ip, T press_Ipp, T vel_dir_Im, T vel_dir_Ip, const bool &gf_is_rho,
    const reconstruct_params_t &reconstruct_params) {

  using std::fabs;
  using std::max;
  using std::min;

  // Unpack all PPM parameters
  const bool &ppm_shock_detection = reconstruct_params.ppm_shock_detection;
  const bool &ppm_zone_flattening = reconstruct_params.ppm_zone_flattening;
  const T &poly_k = reconstruct_params.poly_k;
  const T &poly_gamma = reconstruct_params.poly_gamma;
  const T &ppm_eta1 = reconstruct_params.ppm_eta1;
  const T &ppm_eta2 = reconstruct_params.ppm_eta2;
  const T &ppm_eps = reconstruct_params.ppm_eps;
  const T &ppm_eps_shock = reconstruct_params.ppm_eps_shock;
  const T &ppm_small = reconstruct_params.ppm_small;
  const T &ppm_omega1 = reconstruct_params.ppm_omega1;
  const T &ppm_omega2 = reconstruct_params.ppm_omega2;

  // Helpers
  const T diff_Im = gf_I - gf_Imm;
  const T diff_I = gf_Ip - gf_Im;
  const T diff_Ip = gf_Ipp - gf_I;

  const T delta_Im = 0.5 * diff_Im;
  const T delta_I = 0.5 * diff_I;
  const T delta_Ip = 0.5 * diff_Ip;

  const T _2fabs_Im_Imm = 2 * fabs(gf_Im - gf_Imm);
  const T _2fabs_I_Im = 2 * fabs(gf_I - gf_Im);
  const T _2fabs_Ip_I = 2 * fabs(gf_Ip - gf_I);
  const T _2fabs_Ipp_Ip = 2 * fabs(gf_Ipp - gf_Ip);

  const bool same_sgn_Im = ((gf_I - gf_Im) * (gf_Im - gf_Imm) > 0);
  const bool same_sgn_I = ((gf_Ip - gf_I) * (gf_I - gf_Im) > 0);
  const bool same_sgn_Ip = ((gf_Ipp - gf_Ip) * (gf_Ip - gf_I) > 0);

  const T deltamod_Im =
      same_sgn_Im
          ? sgn(delta_Im) * min(fabs(delta_Im), min(_2fabs_Im_Imm, _2fabs_I_Im))
          : 0;
  const T deltamod_I =
      same_sgn_I
          ? sgn(delta_I) * min(fabs(delta_I), min(_2fabs_I_Im, _2fabs_Ip_I))
          : 0;
  const T deltamod_Ip =
      same_sgn_Ip
          ? sgn(delta_Ip) * min(fabs(delta_Ip), min(_2fabs_Ip_I, _2fabs_Ipp_Ip))
          : 0;

  /* Initial reconstructed states at the interfaces between cells Im/I ans I/Ip
   * NOTE: not const because they may change later */
  const T gf_Im_I =
      0.5 * (gf_Im + gf_I) + (1. / 6.) * (deltamod_Im - deltamod_I);
  const T gf_I_Ip =
      0.5 * (gf_I + gf_Ip) + (1. / 6.) * (deltamod_I - deltamod_Ip);

  T rc_low = gf_Im_I;
  T rc_up = gf_I_Ip;

  const T diff_press_I = press_Ip - press_Im;
  const T min_press_I = min(press_Im, press_Ip);

  /* Shock detection (eqs. 1.15, 1.16, 1.17 with uniform grid spacing).
   * This is only applied to rho and only if the shock is marked as a contact
   * discontinuity (see eq. 3.2) */
  // FIXME: contact discontinuity check only valid for polytropic/ideal-fluid
  // EOS
  if (ppm_shock_detection and gf_is_rho) {
    const T k_gamma = poly_k * poly_gamma;
    const bool contact_disc = (k_gamma * fabs(diff_I) * min_press_I >=
                               fabs(diff_press_I) * min(gf_Ip, gf_Im));
    if (contact_disc) {
      // This assumes gf_var is rho (gf_is_rho is true)
      const T d2rho_Im = gf_I - 2 * gf_Im + gf_Imm;
      const T d2rho_Ip = gf_Ipp - 2 * gf_Ip + gf_I;
      const T d2_prod = d2rho_Im * d2rho_Ip;
      const bool cond2 =
          (fabs(diff_I) - ppm_eps_shock * min(fabs(gf_Ip), fabs(gf_Im))) > 0.;

      const T eta_tilde_I = ((d2_prod < 0) and cond2)
                                ? ((-1. / 6.) * (d2rho_Ip - d2rho_Im) / diff_I)
                                : 0;
      const T eta_I = max(0., min(ppm_eta1 * (eta_tilde_I - ppm_eta2), 1.));

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
    const T w_I = ((fabs(diff_press_I) > ppm_eps * min_press_I) and
                   (vel_dir_Im > vel_dir_Ip))
                      ? 1.
                      : 0;
    const T ftilde_I =
        (fabs(press_Ipp - press_Imm) < ppm_small)
            ? 1.0
            : max(0.,
                  1. - w_I * max(0., ppm_omega2 * ((diff_press_I /
                                                    (press_Ipp - press_Imm)) -
                                                   ppm_omega1)));

    const T one_minus_ftilde_I_gfI = (1 - ftilde_I) * gf_I;

    rc_low = ftilde_I * rc_low + one_minus_ftilde_I_gfI;
    rc_up = ftilde_I * rc_up + one_minus_ftilde_I_gfI;
  }

  // Monotonization (see eq. 1.10)
  if ((rc_up - gf_I) * (gf_I - rc_low) <= 0) {
    rc_low = gf_I;
    rc_up = gf_I;
  } else {
    const T diff_rc = rc_up - rc_low;
    const T diff_rc_sq = diff_rc * diff_rc;
    const T gf6_I = 6 * diff_rc * (gf_I - 0.5 * (rc_low + rc_up));

    if (gf6_I > diff_rc_sq) {
      rc_low = 3 * gf_I - 2 * rc_up;
    }

    else if (gf6_I < -diff_rc_sq) {
      rc_up = 3 * gf_I - 2 * rc_low;
    }
  }

  // Return the lower and upper reconstructed states in cell I
  const array<T, 2> rc = {rc_low, rc_up};
  return rc;
}

template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<T, 2>
ppm_reconstruct(T gf_Immm, T gf_Imm, T gf_Im, T gf_Ip, T gf_Ipp, T gf_Ippp,
                T press_Immm, T press_Imm, T press_Im, T press_Ip, T press_Ipp,
                T press_Ippp, T vel_dir_Imm, T vel_dir_Im, T vel_dir_Ip,
                T vel_dir_Ipp, const bool &gf_is_rho,
                const reconstruct_params_t &reconstruct_params) {

  const auto rc_Im{ppm(gf_Immm, gf_Imm, gf_Im, gf_Ip, gf_Ipp, press_Immm,
                       press_Imm, press_Ip, press_Ipp, vel_dir_Imm, vel_dir_Ip,
                       gf_is_rho, reconstruct_params)};

  const auto rc_Ip{ppm(gf_Imm, gf_Im, gf_Ip, gf_Ipp, gf_Ippp, press_Imm,
                       press_Im, press_Ipp, press_Ippp, vel_dir_Im, vel_dir_Ipp,
                       gf_is_rho, reconstruct_params)};

  return {rc_Im[1], rc_Ip[0]};
}

} // namespace ReconX

#endif // RECONX_PPM_HXX
