#ifndef C2P_1DREPRIMAND_INTERNALS_HXX
#define C2P_1DREPRIMAND_INTERNALS_HXX

#include <cmath>
#include <limits>

#include "prims.hxx"
#include "cons.hxx"
#include "c2p_1DRePrimAnd_intervals.hxx"
#include "c2p_1DRePrimAnd_rootfinder.hxx"

namespace Con2PrimFactory {
namespace RePrimAnd {

// --------------------------- f_upper ---------------------------------

class f_upper {
public:
  using value_t = CCTK_REAL;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE
  f_upper(CCTK_REAL h0_, CCTK_REAL rsqr_, CCTK_REAL rbsqr_, CCTK_REAL bsqr_)
      : h0(h0_), h0sqr(h0_ * h0_), rsqr(rsqr_), rbsqr(rbsqr_), bsqr(bsqr_) {}

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REAL
  x_from_mu(CCTK_REAL mu) const {
    return 1.0 / (1.0 + mu * bsqr);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REAL
  rfsqr_from_mu_x(CCTK_REAL mu, CCTK_REAL x) const {
    return x * (rsqr * x + mu * (x + 1.0) * rbsqr);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REAL
  new_h0w_from_mu_x(CCTK_REAL mu, CCTK_REAL x) const {
    return std::sqrt(h0sqr + rfsqr_from_mu_x(mu, x));
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE auto
  operator()(CCTK_REAL mu) const -> std::pair<CCTK_REAL, CCTK_REAL> {
    const CCTK_REAL x = x_from_mu(mu);
    const CCTK_REAL xsq = x * x;
    const CCTK_REAL hw = new_h0w_from_mu_x(mu, x);
    const CCTK_REAL b = x * (xsq * rsqr + mu * (1.0 + x + xsq) * rbsqr);
    const CCTK_REAL f = mu * hw - 1.0;
    const CCTK_REAL df = (h0sqr + b) / hw;
    return {f, df};
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE interval<CCTK_REAL>
  initial_bracket() const {
    CCTK_REAL mu_min = 1.0 / std::sqrt(h0sqr + rsqr);
    const CCTK_REAL mu0 = 1.0 / h0;
    const CCTK_REAL rfs_min = rfsqr_from_mu_x(mu0, x_from_mu(mu0));
    CCTK_REAL mu_max = 1.0 / std::sqrt(h0sqr + rfs_min);
    const CCTK_REAL margin = 10 * std::numeric_limits<CCTK_REAL>::epsilon();
    mu_max *= 1.0 + margin;
    mu_min *= 1.0 - margin;
    if (mu_max <= mu_min) {
      mu_min = 0.0;
      mu_max = mu0 * (1.0 + margin);
    }
    return {mu_min, mu_max};
  }

  CCTK_REAL h0, h0sqr, rsqr, rbsqr, bsqr;
};

// --------------------------- froot -----------------------------------

template <typename EOSType> class froot {
public:
  using value_t = CCTK_REAL;

  struct cache {
    CCTK_REAL ye{0}, lmu{0}, x{0};
    CCTK_REAL rho{0}, rho_raw{0};
    CCTK_REAL eps{0}, eps_raw{0};
    CCTK_REAL etot{0}, etot_raw{0};
    CCTK_REAL press{0};
    CCTK_REAL vsqr{0};
    CCTK_REAL w{1};
    unsigned int calls{0};
  };

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE
  froot(const EOSType *eos_, CCTK_REAL valid_ye, CCTK_REAL d_, CCTK_REAL qtot_,
        CCTK_REAL rsqr_, CCTK_REAL rbsqr_, CCTK_REAL bsqr_, cache &last_)
      : eos(eos_), d(d_), qtot(qtot_), rsqr(rsqr_), rbsqr(rbsqr_), bsqr(bsqr_),
        brosqr(rsqr_ * bsqr_ - rbsqr_), last(last_) {
    last.ye = valid_ye;
    last.calls = 0;

    // For P >= 0, h >= 1 + eps_min over the entire EOS domain. Sampling
    // h at rho_min is not a global bound; negative physical eps permits h < 1.
    // solve checks eps_min > -1 before constructing the root function.
    h0 = 1.0 + fmin(CCTK_REAL(0.0), eos->rgeps.min);

    const CCTK_REAL zsqrinf = rsqr / (h0 * h0);
    const CCTK_REAL wsqrinf = 1.0 + zsqrinf;
    winf = std::sqrt(wsqrinf);
    vsqrinf = zsqrinf / wsqrinf;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REAL
  operator()(CCTK_REAL mu) {
    cache &c = last;
    c.lmu = mu;
    c.x = x_from_mu(mu);

    const CCTK_REAL rfsqr = rfsqr_from_mu_x(mu, c.x);
    const CCTK_REAL qf = qf_from_mu_x(mu, c.x);

    c.vsqr = rfsqr * mu * mu;
    if (c.vsqr >= vsqrinf) {
      c.vsqr = vsqrinf;
      c.w = winf;
    } else {
      c.w = 1.0 / std::sqrt(1.0 - c.vsqr);
    }

    c.rho_raw = d / c.w;
    c.rho = fmin(fmax(eos->rgrho.min, c.rho_raw), eos->rgrho.max);

    c.eps_raw = get_eps_raw(mu, qf, rfsqr, c.w);
    const auto rgeps = eos->range_eps_from_rho_ye(c.rho, c.ye);
    c.eps = fmin(fmax(rgeps.min, c.eps_raw), rgeps.max);

    c.etot_raw = c.rho_raw * (1.0 + c.eps_raw);
    c.etot = c.rho * (1.0 + c.eps);

    c.press = eos->press_from_rho_eps_ye(c.rho, c.eps, c.ye);
    ++c.calls;

    const CCTK_REAL inv_etot =
        CCTK_REAL(1.0) /
        fmax(c.etot, std::numeric_limits<CCTK_REAL>::min());
    const CCTK_REAL one_plus_a = (c.etot + c.press) * inv_etot;
    const CCTK_REAL h = (c.etot + c.press) / c.rho;

    const CCTK_REAL hbw_raw = one_plus_a * (1.0 + qf - mu * rfsqr);
    const CCTK_REAL hbw = fmax(hbw_raw, h / c.w);
    const CCTK_REAL newmu = 1.0 / (hbw + rfsqr * mu);

    return mu - newmu;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE interval<CCTK_REAL>
  initial_bracket(ROOTSTAT &status) const {
    status = ROOTSTAT::SUCCESS;
    CCTK_REAL mu_max = CCTK_REAL(1) / h0;

    if (rsqr >= h0 * h0) {
      const int ndigits2 = 36;
      const CCTK_REAL margin = std::ldexp(CCTK_REAL(1), 3 - ndigits2);
      RePrimAnd::f_upper g(h0, rsqr, rbsqr, bsqr);
      const CCTK_REAL mu_try =
          findroot_using_deriv(g, status, ndigits2, ndigits2 + 4);

      if (status != ROOTSTAT::SUCCESS || !std::isfinite(mu_try)) {
        if (status == ROOTSTAT::SUCCESS)
          status = ROOTSTAT::NOT_CONVERGED;
        return {CCTK_REAL(0), CCTK_REAL(1) / h0};
      }

      mu_max = mu_try * (CCTK_REAL(1) + margin);
    }

    return {CCTK_REAL(0), mu_max};
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE interval<CCTK_REAL>
  initial_bracket() const {
    ROOTSTAT status{};
    return initial_bracket(status);
  }

  // RePrimAnd stop condition used by the TOMS748-style bracket solver.
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE bool
  stopif(CCTK_REAL mu, CCTK_REAL dmu, CCTK_REAL acc) const {
    return std::abs(dmu) * last.w * last.w < mu * acc;
  }

  // Public parameters accessed by rarecase/f_rare
  CCTK_REAL h0{1}, d{0}, qtot{0}, rsqr{0}, rbsqr{0}, bsqr{0}, brosqr{0};
  CCTK_REAL winf{1}, vsqrinf{0};

private:
  const EOSType *eos;
  cache &last;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REAL
  x_from_mu(CCTK_REAL mu) const {
    return 1.0 / (1.0 + mu * bsqr);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REAL
  rfsqr_from_mu_x(CCTK_REAL mu, CCTK_REAL x) const {
    return x * (rsqr * x + mu * (x + 1.0) * rbsqr);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REAL
  qf_from_mu_x(CCTK_REAL mu, CCTK_REAL x) const {
    const CCTK_REAL mux = mu * x;
    return qtot - (bsqr + mux * mux * brosqr) * CCTK_REAL(0.5);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE static CCTK_REAL
  get_eps_raw(CCTK_REAL mu, CCTK_REAL qf, CCTK_REAL rfsqr, CCTK_REAL w) {
    return w * (qf - mu * rfsqr * (1.0 - mu * w / (1.0 + w)));
  }
};

// -------------------- f_rare and rarecase -------------------

template <typename EOSType> class f_rare {
  const CCTK_REAL v2targ;
  const froot<EOSType> &f;

public:
  using value_t = CCTK_REAL;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE
  f_rare(CCTK_REAL wtarg, const froot<EOSType> &f_)
      : v2targ(1.0 - 1.0 / (wtarg * wtarg)), f(f_) {}

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE auto
  operator()(CCTK_REAL mu) const -> std::pair<CCTK_REAL, CCTK_REAL> {
    const CCTK_REAL x = 1.0 / (1.0 + mu * f.bsqr);
    const CCTK_REAL xsq = x * x;
    const CCTK_REAL rfsq = x * (f.rsqr * x + mu * (x + 1.0) * f.rbsqr);
    const CCTK_REAL vsq = mu * mu * rfsq;
    const CCTK_REAL y = vsq - v2targ;
    const CCTK_REAL dy =
        CCTK_REAL(2) * mu * x * (xsq * f.rsqr + mu * (xsq + x + 1.0) * f.rbsqr);
    return {y, dy};
  }
};

template <typename EOSType> class rarecase {
public:
  interval<CCTK_REAL> bracket;
  bool rho_too_big{false}, rho_big{false}, rho_too_small{false},
      rho_small{false};

  CCTK_HOST CCTK_DEVICE rarecase(const interval<CCTK_REAL> ibracket,
                                 const interval<CCTK_REAL> rgrho,
                                 const froot<EOSType> &f) {
    CCTK_REAL muc0 = ibracket.min();
    CCTK_REAL muc1 = ibracket.max();
    const int ndigits = 30;

    if (f.d > rgrho.max()) {
      const CCTK_REAL wc = f.d / rgrho.max();
      if (wc > f.winf) {
        rho_too_big = true;
      } else {
        f_rare<EOSType> g(wc, f);
        if (g(muc1).first <= 0.0) {
          rho_too_big = true;
        } else if (g(muc0).first < 0.0) {
          ROOTSTAT status{};
          CCTK_REAL mucc =
              findroot_using_deriv(g, ibracket, status, ndigits, ndigits + 2);
          if (status == ROOTSTAT::SUCCESS) {
            muc0 = fmax(muc0, mucc);
            rho_big = true;
          } else {
            rho_too_big = true;
          }
        }
      }
    }

    if (f.d < f.winf * rgrho.min()) {
      const CCTK_REAL wc = f.d / rgrho.min();
      if (wc < 1.0) {
        rho_too_small = true;
      } else {
        f_rare<EOSType> g(wc, f);
        if (g(muc0).first >= 0.0) {
          rho_too_small = true;
        } else if (g(muc1).first > 0.0) {
          ROOTSTAT status{};
          CCTK_REAL mucc =
              findroot_using_deriv(g, ibracket, status, ndigits, ndigits + 2);
          if (status == ROOTSTAT::SUCCESS) {
            muc1 = fmin(muc1, mucc);
            rho_small = true;
          } else {
            rho_too_small = true;
          }
        }
      }
    }

    bracket = interval<CCTK_REAL>{muc0, muc1};
  }
};

} // namespace RePrimAnd
} // namespace Con2PrimFactory

#endif
