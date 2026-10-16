#ifndef EOS_3P_TABULATED3D_HXX
#define EOS_3P_TABULATED3D_HXX

#include <cmath>
#include <cassert>
#include <limits>
#include <string>
#include <array>
#include <mpi.h>
#include <hdf5.h>

#include "../eos_3p.hxx"
#include "../utils/eos_brent.hxx" // zero_brent
#include "../utils/eos_linear_interp_ND.hxx"

#ifndef NTABLES
#define NTABLES 19
#endif

namespace EOSX {
using namespace std;
using namespace eos_constants;

//------------------------------------
// Reader-agnostic raw table container
//------------------------------------
struct eos_tabulated3d_raw_table_t {
  int nrho = 0;
  int ntemp = 0;
  int nye = 0;
  int npoints = 0;
  int have_rel_cs2 = 0;

  // Final device/managed storage
  CCTK_REAL *logrho = nullptr;
  CCTK_REAL *logtemp = nullptr;
  CCTK_REAL *yes = nullptr;
  CCTK_REAL *alltables = nullptr;
  CCTK_REAL *energy_shift = nullptr;
};

CCTK_HOST void eos_readtable(const std::string &filename,
                             eos_tabulated3d_raw_table_t &tab);

//-----------------
// Tabulated 3D EOS
//-----------------
class eos_3p_tabulated3d : public eos_3p {
public:
  // must match order of HDF5 datasets below
  enum EV {
    PRESS = 0,   // "logpress"
    EPS = 1,     // "logenergy"
    S = 2,       // "entropy"
    MUNU = 3,    // "munu"
    CS2 = 4,     // "cs2"
    DEDT = 5,    // "dedt"
    DPDRHOE = 6, // "dpdrhoe"
    DPDERHO = 7, // "dpderho"
    MUHAT = 8,   // "muhat"
    MU_E = 9,    // "mu_e"
    MU_P = 10,   // "mu_p"
    MU_N = 11,   // "mu_n"
    XA = 12,     // "Xa"
    XH = 13,     // "Xh"
    XN = 14,     // "Xn"
    XP = 15,     // "Xp"
    ABAR = 16,   // "Abar"
    ZBAR = 17,   // "Zbar"
    GAMMA = 18,  // "Gamma"
    NUM_VARS
  };

  CCTK_REAL gamma; // table Γ
  CCTK_REAL *energy_shift;
  range rgeps;
  linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> *interptable;

  CCTK_HOST void init(const std::string &filename, range &rgeps_out,
                      const range &, const range &) {

    // Read raw table
    eos_tabulated3d_raw_table_t tab;
    eos_readtable(filename, tab);

    const int nrho = tab.nrho;
    const int ntemp = tab.ntemp;
    const int nye = tab.nye;
    const int npoints = tab.npoints;

    const int have_rel_cs2 = tab.have_rel_cs2;

    CCTK_REAL *logrho = tab.logrho;
    CCTK_REAL *logtemp = tab.logtemp;
    CCTK_REAL *yes = tab.yes;
    CCTK_REAL *alltables = tab.alltables;
    energy_shift = tab.energy_shift;

    // -----------------------------------------------------------------
    // Common post-processing (unit conversions, log conversions, interp)
    // -----------------------------------------------------------------

    *energy_shift *= EPSGF;
    const CCTK_REAL ln10 = log(10.0);
    const CCTK_REAL inv_time2 = 1 / (TIMEGF * TIMEGF);
    for (int i = 0; i < nrho; i++)
      logrho[i] = logrho[i] * ln10 + log(RHOGF);
    for (int i = 0; i < ntemp; i++)
      logtemp[i] *= ln10;

    // cs2 handling:
    // - Convert cs2 to code units
    // - If NOT already relativistic, divide by h
    // - Clamp to <= 0.9999999
    const CCTK_REAL max_cs2 = 0.9999999;

    for (size_t idx = 0; idx < (size_t)npoints; idx++) {
      size_t b = idx * NTABLES;

      alltables[b + PRESS] = alltables[b + PRESS] * ln10 + log(PRESSGF);
      alltables[b + EPS] = alltables[b + EPS] * ln10 + log(EPSGF);

      // Convert cs2 to code units first
      CCTK_REAL cs2 = alltables[b + CS2] * LENGTHGF * LENGTHGF * inv_time2;
      if (cs2 < 0)
        cs2 = 0;

      // If cs2 is not already relativistic, divide by enthalpy h
      if (!have_rel_cs2) {
        const int irho = (int)(idx % (size_t)nrho);
        const CCTK_REAL rhoL = exp(logrho[irho]);

        const CCTK_REAL pressL = exp(alltables[b + PRESS]);

        const CCTK_REAL epspL = exp(alltables[b + EPS]); // eps + shift
        const CCTK_REAL epsL = epspL - *energy_shift;

        const CCTK_REAL hL = 1.0 + epsL + pressL / rhoL;
        cs2 /= hL;
      }

      if (cs2 > max_cs2)
        cs2 = max_cs2;
      alltables[b + CS2] = cs2;

      alltables[b + DEDT] = alltables[b + DEDT] * EPSGF;
      alltables[b + DPDRHOE] = alltables[b + DPDRHOE] * PRESSGF / RHOGF;
      alltables[b + DPDERHO] = alltables[b + DPDERHO] * PRESSGF / EPSGF;
    }

    std::array<size_t, 3> dims{size_t(nrho), size_t(ntemp), size_t(nye)};
    interptable = (decltype(interptable))amrex::The_Managed_Arena()->alloc(
        sizeof(*interptable));
    new (interptable) linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES>(
        alltables, dims, logrho, logtemp, yes);

    set_range_rho(
        range{exp(interptable->xmin<0>()), exp(interptable->xmax<0>())});
    set_range_temp(
        range{exp(interptable->xmin<1>()), exp(interptable->xmax<1>())});
    set_range_ye(range{interptable->xmin<2>(), interptable->xmax<2>()});

    rgeps = compute_eps_range_full_table();
    rgeps_out = rgeps;
  }

  // Device-callable routines

  //// Table inversions
  template <size_t var>
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  logtemp_from_rho_var_ye(const CCTK_REAL rho, CCTK_REAL &invar,
                          const CCTK_REAL ye) const {
    // bound inputs
    CCTK_REAL r = std::fmin(std::fmax(rho, rgrho.min), rgrho.max);
    CCTK_REAL lrho = std::log(r);

    // table-edge clamp
    auto vmin =
        interptable->interpolate<var>(lrho, interptable->xmin<1>(), ye)[0];
    auto vmax =
        interptable->interpolate<var>(lrho, interptable->xmax<1>(), ye)[0];
    if (invar <= vmin) {
      invar = vmin;
      return interptable->xmin<1>();
    }
    if (invar >= vmax) {
      invar = vmax;
      return interptable->xmax<1>();
    }

    // root-find for logtemp
    auto func = [&](CCTK_REAL &lt) {
      CCTK_REAL val = interptable->interpolate<var>(lrho, lt, ye)[0];
      return invar - val;
    };
    return zero_brent(interptable->xmin<1>(), interptable->xmax<1>(), 1.e-14,
                      func);
  }

  template<size_t var>
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  logrho_from_var_temp_ye(CCTK_REAL &invar, const CCTK_REAL temp, 
                                const CCTK_REAL ye) const {
    // bound inputs
    CCTK_REAL t = std::fmin(std::fmax(temp, rgtemp.min), rgtemp.max);
    CCTK_REAL lt = std::log(t);

    // table‐edge clamp
    auto vmin =
        interptable->interpolate<var>(interptable->xmin<0>(), lt, ye)[0];
    auto vmax =
        interptable->interpolate<var>(interptable->xmax<0>(), lt, ye)[0];
    if (invar <= vmin) {
      invar = vmin;
      return interptable->xmin<0>();
    }
    if (invar >= vmax) {
      invar = vmax;
      return interptable->xmax<0>();
    }

    // root‐find for logrho
    auto func = [&](CCTK_REAL &lrho) {
      CCTK_REAL val = interptable->interpolate<var>(lrho, lt, ye)[0];
      return invar - val;
    };
    return zero_brent(interptable->xmin<0>(), interptable->xmax<0>(), 1.e-14,
                      func);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  logtemp_from_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                          const CCTK_REAL ye) const {
    // Bound physical eps before taking the logarithm. The inverse below
    // applies the local temperature-edge bounds at this (rho, Ye).
    eps = std::fmin(std::fmax(eps, rgeps.min), rgeps.max);
    const CCTK_REAL shifted_eps = eps + *energy_shift;
    assert(shifted_eps > 0.0);
    CCTK_REAL leps = std::log(shifted_eps);
    CCTK_REAL lt = logtemp_from_rho_var_ye<EV::EPS>(rho, leps, ye);
    eps = exp(leps) - *energy_shift;
    return lt;
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  temp_from_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                       const CCTK_REAL ye) const {
    CCTK_REAL lt = logtemp_from_rho_eps_ye(rho, eps, ye);
    return exp(lt);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  temp_from_rho_entropy_ye(const CCTK_REAL rho, CCTK_REAL &ent,
                           const CCTK_REAL ye) const {
    CCTK_REAL lt = logtemp_from_rho_var_ye<EV::S>(rho, ent, ye);
    return exp(lt);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  temp_from_rho_press_ye(const CCTK_REAL rho, CCTK_REAL &press,
                                const CCTK_REAL ye) const {
    CCTK_REAL lP = log(press);
    CCTK_REAL lt = logtemp_from_rho_var_ye<EV::PRESS>(rho, lP, ye);
    press = exp(lP);
    return exp(lt);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  rho_from_press_temp_ye(CCTK_REAL &press, const CCTK_REAL temp,
                                const CCTK_REAL ye) const {
    CCTK_REAL lP = log(press);
    CCTK_REAL lr = logrho_from_var_temp_ye<EV::PRESS>(lP, temp, ye);
    press = exp(lP);
    return exp(lr);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  press_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                         const CCTK_REAL ye) const {
    // bound
    CCTK_REAL r = std::fmin(std::fmax(rho, rgrho.min), rgrho.max);
    CCTK_REAL t = std::fmin(std::fmax(temp, rgtemp.min), rgtemp.max);
    CCTK_REAL lr = std::log(r), lt = std::log(t);
    CCTK_REAL v = interptable->interpolate<EV::PRESS>(lr, lt, ye)[0];
    return exp(v);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  press_from_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                        const CCTK_REAL ye) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL lt = logtemp_from_rho_eps_ye(rho, eps, ye);
    CCTK_REAL v = interptable->interpolate<EV::PRESS>(lr, lt, ye)[0];
    return exp(v);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  eps_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                       const CCTK_REAL ye) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL lt = std::log(std::fmin(std::fmax(temp, rgtemp.min), rgtemp.max));
    CCTK_REAL v = interptable->interpolate<EV::EPS>(lr, lt, ye)[0];
    return exp(v) - *energy_shift;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_rho_press_ye(const CCTK_REAL rho, const CCTK_REAL press,
                        const CCTK_REAL ye) const {

    printf(
        "This routine should not be used. There is no monotonicity condition "
        "to enforce a succesfull inversion from eps(press). So you better "
        "rewrite your code to not require this call. \n");
    assert(false);
    return 0;
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  csnd_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                        const CCTK_REAL ye) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL lt = std::log(std::fmin(std::fmax(temp, rgtemp.min), rgtemp.max));
    CCTK_REAL v = interptable->interpolate<EV::CS2>(lr, lt, ye)[0];
    if (v < 0) {
      printf("cs^2 < 0 detected! This should have been fixed by table "
             "preprocessing!\n");
      v = 0;
    }
    return sqrt(v);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  csnd_from_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                       const CCTK_REAL ye) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL lt = logtemp_from_rho_eps_ye(rho, eps, ye);
    CCTK_REAL v = interptable->interpolate<EV::CS2>(lr, lt, ye)[0];
    if (v < 0) {
      printf("cs^2 < 0 detected! This should have been fixed by table "
             "preprocessing!\n");
      v = 0;
    }
    return sqrt(v);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  press_derivs_from_rho_eps_ye(CCTK_REAL &press, CCTK_REAL &dpdrho,
                               CCTK_REAL &dpdeps, const CCTK_REAL rho,
                               const CCTK_REAL eps, const CCTK_REAL ye) const {
    printf("press_derivs_from_rho_eps_ye is not supported for now! \n");
    assert(false);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  entropy_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                           const CCTK_REAL ye) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL lt = std::log(std::fmin(std::fmax(temp, rgtemp.min), rgtemp.max));
    return interptable->interpolate<EV::S>(lr, lt, ye)[0];
  }

  CCTK_HOST CCTK_DEVICE inline void
  mu_pne_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                          const CCTK_REAL ye, CCTK_REAL &mup, CCTK_REAL &mun,
                          CCTK_REAL &mue) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL lt = std::log(std::fmin(std::fmax(temp, rgtemp.min), rgtemp.max));
    mup = interptable->interpolate<EV::MU_P>(lr, lt, ye)[0];
    mun = interptable->interpolate<EV::MU_N>(lr, lt, ye)[0];
    mue = interptable->interpolate<EV::MU_E>(lr, lt, ye)[0];
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  mu_lepton_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                             const CCTK_REAL ye) const {
    CCTK_REAL mup, mun, mue;
    mu_pne_from_rho_temp_ye(rho, temp, ye, mup, mun, mue);
    return mue + mup - mun;
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  ye_beq_from_rho_temp(const CCTK_REAL rho, const CCTK_REAL temp) const {
    // Neutrino-less beta equilibrium: mu_e + mu_p - mu_n = 0.
    const CCTK_REAL lr =
        log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    const CCTK_REAL lt =
        log(std::fmin(std::fmax(temp, rgtemp.min), rgtemp.max));
    const auto func = [&](const CCTK_REAL Ye) {
      const auto mu = interptable->interpolate<EV::MU_E, EV::MU_P, EV::MU_N>(
          lr, lt, Ye);
      return mu[0] + mu[1] - mu[2];
    };

    const auto *yes = interptable->x[2];
    size_t a = 0, b = interptable->num_points[2] - 1;
    CCTK_REAL fa = func(yes[a]), fb = func(yes[b]);
    assert(std::isfinite(fa) && std::isfinite(fb));
    if (fa == 0.0)
      return yes[a];
    if (fb == 0.0)
      return yes[b];

    // If equilibrium lies outside the table, use the closest endpoint.
    if ((fa < 0.0) == (fb < 0.0))
      return fabs(fa) <= fabs(fb) ? yes[a] : yes[b];

    // Locate the sign-changing Ye cell. The chemical potentials are linear
    // in Ye inside this cell, so the final interpolation gives its root.
    while (b - a > 1) {
      const size_t m = a + (b - a) / 2;
      const CCTK_REAL fm = func(yes[m]);
      assert(std::isfinite(fm));
      if (fm == 0.0)
        return yes[m];
      if ((fa < 0.0) != (fm < 0.0)) {
        b = m;
        fb = fm;
      } else {
        a = m;
        fa = fm;
      }
    }

    const CCTK_REAL scale = std::fmax(fabs(fa), fabs(fb));
    const CCTK_REAL wa = fabs(fa) / scale, wb = fabs(fb) / scale;
    return yes[a] + (yes[b] - yes[a]) * (wa / (wa + wb));
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_rho_kappa_ye(const CCTK_REAL rho,
                          const CCTK_REAL kappa, // kappa=entropy
                          const CCTK_REAL ye) const {
    printf("press_from_rho_kappa_ye is not supported for tabulated EOS! \n");
    assert(false);
    return 0.0;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_rho_kappa_ye(const CCTK_REAL rho,
                        const CCTK_REAL kappa, // kappa=entropy
                        const CCTK_REAL ye) const {
    printf("eps_from_rho_kappa_ye is not supported for tabulated EOS! \n");
    assert(false);
    return 0.0;
  };

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  kappa_from_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                        const CCTK_REAL ye) const {
    return entropy_from_rho_temp_ye(rho, temp_from_rho_eps_ye(rho, eps, ye),
                                    ye);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  kappa_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                        const CCTK_REAL ye) const {
    // Table entropy is evolved kappa. Reuse the known temperature.
    return entropy_from_rho_temp_ye(rho, temp, ye);
  }

  CCTK_HOST CCTK_DEVICE inline range
  range_eps_from_rho_ye(const CCTK_REAL rho, const CCTK_REAL ye) const {
    CCTK_REAL lr = std::log(std::fmin(std::fmax(rho, rgrho.min), rgrho.max));
    CCTK_REAL vmin =
        interptable->interpolate<EV::EPS>(lr, interptable->xmin<1>(), ye)[0];
    CCTK_REAL vmax =
        interptable->interpolate<EV::EPS>(lr, interptable->xmax<1>(), ye)[0];
    return range{exp(vmin) - *energy_shift, exp(vmax) - *energy_shift};
  }

  CCTK_HOST CCTK_DEVICE inline range compute_eps_range_full_table() const {
    size_t n0 = interptable->num_points[0];
    size_t n1 = interptable->num_points[1];
    size_t n2 = interptable->num_points[2];
    size_t total = n0 * n1 * n2;

    CCTK_REAL eps_min = std::numeric_limits<CCTK_REAL>::max();
    CCTK_REAL eps_max = std::numeric_limits<CCTK_REAL>::lowest();

    for (size_t i = 0; i < total; i++) {
      const CCTK_REAL logeps = interptable->y[EV::EPS + NTABLES * i];
      // Only interpolation uses shifted energy; callers use physical eps.
      CCTK_REAL val = exp(logeps) - *energy_shift;
      eps_min = std::fmin(eps_min, val);
      eps_max = std::fmax(eps_max, val);
    }

    return range{eps_min, eps_max};
  }

}; // class eos_3p_tabulated3d

} // namespace EOSX

#endif // EOS_3P_TABULATED3D_HXX
