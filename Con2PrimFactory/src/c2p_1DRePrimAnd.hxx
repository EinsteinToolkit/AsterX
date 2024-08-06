/*! \file c2p_1DRePrimAnd.hxx 
\brief Ideal MHD primitive variable recovery algorithm
*/

#ifndef C2P_1DREPRIMAND_HXX
#define C2P_1DREPRIMAND_HXX

#include <cassert>
#include <cmath>
#include <limits>

#include "c2p_1DRePrimAnd_rootfinder.hxx"
#include "c2p_1DRePrimAnd_internals.hxx"

class c2p_1DRePrimAnd : public c2p {

public:

using namespace std;

using range   = eos::range;
using report  = c2p_report;






con2prim_mhd::con2prim_mhd(eos_thermal eos_, CCTK_REAL rho_strict_, 
    bool ye_lenient_, CCTK_REAL z_lim_, CCTK_REAL b_lim_, 
    const atmosphere& atmo_, CCTK_REAL acc_, int max_iter_) 
: eos(std::move(eos_)), rho_strict(rho_strict_), 
  ye_lenient(ye_lenient_), z_lim(z_lim_),
  bsqr_lim(b_lim_*b_lim_), atmo(atmo_), acc(acc_), max_iter(max_iter_)
{
  w_lim = sqrt(1.0 + z_lim*z_lim);
  v_lim = z_lim / w_lim;
}

/**
Requires basic physical constraints
\f[ 
h_0>0, \qquad
r^2 \ge 0, \qquad
(r^l b_l)^2 \ge 0, \qquad
b^2 \ge 0
\f]
The safety margin for the root bracketing needs to satisfy
\f$ 1_m > 1 \f$.
**/
f_upper::f_upper(CCTK_REAL h0_, CCTK_REAL rsqr_, 
  CCTK_REAL rbsqr_, CCTK_REAL bsqr_)
: h0(h0_), h0sqr(h0_*h0_), rsqr(rsqr_), rbsqr(rbsqr_), bsqr(bsqr_)
{
  assert(h0 > 0);
  assert(rsqr >= 0);
  assert(rbsqr >= 0);
  assert(bsqr >= 0);  
}

/**
Computes 
\f[ x = \frac{y}{y+b^2} = \frac{1}{1 + \mu b^2} \f]
**/
CCTK_REAL f_upper::x_from_mu(const CCTK_REAL mu) const
{
  return 1 / (1 + mu * bsqr);
}

/**
Uses the formula
\f[
\bar{r}^2 = x^2 r^2_\perp + r^2_\parallel 
= x \left( r^2 x + \mu \left( x + 1 \right) \left(r^l b_l\right)^2 \right)
\f]
**/
CCTK_REAL f_upper::rfsqr_from_mu_x(const CCTK_REAL mu,
  const CCTK_REAL x) const
{
  return x * (rsqr * x + mu * (x + 1.0) * rbsqr);  
}

/**
Based on the identity 
\f[
hW = \sqrt{h^2 + \bar{r}^2}
\f]
and assuming the minimal value for the enthalpy \f$ h = h_0 \f$
**/
CCTK_REAL f_upper::new_h0w_from_mu_x(const CCTK_REAL mu,
  const CCTK_REAL x) const
{
  return sqrt(h0sqr + rfsqr_from_mu_x(mu,x));  
}


/**
This implements the auxiliary root function as defined in the article
https://doi.org/10.1103/PhysRevD.103.023018
**/
auto f_upper::operator()(const CCTK_REAL mu) const 
-> std::pair<CCTK_REAL, CCTK_REAL>
{
  CCTK_REAL x     = x_from_mu(mu);
  CCTK_REAL xsqr  = x*x;
  CCTK_REAL hw    = new_h0w_from_mu_x(mu, x);
  CCTK_REAL b     = x * (xsqr * rsqr + mu * (1 + x + xsqr) * rbsqr);
  CCTK_REAL f     = mu * hw - 1.;
  CCTK_REAL df    = (h0sqr + b) / hw;
  return {f, df};
}


/**
This computes an initial bracket for the auxiliary root
**/
auto f_upper::initial_bracket() const -> interval<CCTK_REAL>
{
  CCTK_REAL mu_min    = 1. / sqrt(h0sqr + rsqr);
  CCTK_REAL mu0       = 1. / h0;
  CCTK_REAL rfsqr_min = rfsqr_from_mu_x(mu0, x_from_mu(mu0));
  CCTK_REAL mu_max    = 1.0 / sqrt(h0sqr + rfsqr_min);
  //mathematically, mu_max >= mu_min as long as r^2,b^2, r^l b_l
  //respect Schwarz inequality. However, it is possible that 
  // mu_min=mu_max, which is why we have to account for roundoff 
  // errors
  CCTK_REAL margin    = 10*std::numeric_limits<CCTK_REAL>::epsilon();
  mu_max *= 1.0 + margin;
  mu_min *= 1.0 - margin;
  // If this ad-hoc margin was not enough, we just use a much wider 
  // bracket (the original one given in the article).
  if (mu_max <= mu_min) 
  { 
    mu_min = 0;
    mu_max = mu0 * (1.0 + margin);
  }
  assert(mu_max > mu_min);
  return {mu_min, mu_max};
}


/**
This sets the parameters and EOS defining the root function. We also 
compute the electron fraction limited to the allowed range of the EOS.
Further, we compute an upper limit for the velocity from
\f[
z = \frac{\bar{r}}{h} 
  \le \frac{r}{h} \le \frac{r}{h_0} 
\f]
where we used \f$ \bar{r} \le r \f$ and the minimum enthalpy \f$ h_0 \f$
provided by the EOS. 
**/
froot::froot(const eos_thermal& eos_, CCTK_REAL valid_ye,
      CCTK_REAL d_, CCTK_REAL qtot_, CCTK_REAL rsqr_, CCTK_REAL rbsqr_,
      CCTK_REAL bsqr_, cache& last_ )
: eos(eos_), h0(eos_.minimal_h()), 
  rho_range(eos_.range_rho()), d(d_), qtot(qtot_), rsqr(rsqr_), 
  rbsqr(rbsqr_), bsqr(bsqr_), 
  brosqr(rsqr_ * bsqr_ - rbsqr_), last(last_)
{
  assert(eos.range_ye().contains(valid_ye));
  last.ye    = valid_ye;
  last.calls = 0;
  
  CCTK_REAL zsqrinf = rsqr / (h0*h0);
  CCTK_REAL wsqrinf = 1 + zsqrinf;
  winf    = sqrt(wsqrinf);
  vsqrinf = zsqrinf / wsqrinf;
}


/**
Computes 
\f[ x = \frac{y}{y+b^2} = \frac{1}{1 + \mu b^2} \f]
**/
CCTK_REAL froot::x_from_mu(const CCTK_REAL mu) const
{
  return 1 / (1 + mu * bsqr);  
}

/**
Uses the formula
\f[
\bar{r}^2 = x^2 r^2_\perp + r^2_\parallel 
= x \left( r^2 x + \mu \left( x + 1 \right) \left(r^l b_l\right)^2 \right)
\f]
**/
CCTK_REAL froot::rfsqr_from_mu_x(const CCTK_REAL mu, 
  const CCTK_REAL x) const
{
  return x * (rsqr * x + mu * (x + 1.0) * rbsqr);  
}

/**
Uses the formula
\f[
\bar{q} = q - \frac{1}{2} \left( b^2 + \mu^2 x^2 b^2 r^2_\perp \right)
\f]
**/
CCTK_REAL froot::qf_from_mu_x(const CCTK_REAL mu, 
  const CCTK_REAL x) const
{
  CCTK_REAL mux = mu * x;
  return qtot - (bsqr + mux*mux*brosqr) / 2;
}
/**
The formula is written in a way that is accurate also for small velocities.
\f[
\epsilon = W \left( \bar{q} - \mu \bar{r}^2 \right) 
           + v^2 \frac{W^2}{1 + W}
\f]
**/
CCTK_REAL froot::get_eps_raw(const CCTK_REAL mu, const CCTK_REAL qf, 
  const CCTK_REAL rfsqr, const CCTK_REAL w) 
{
  return w * (qf - mu * rfsqr*(1.0 - mu * w / (1 + w)));
}

/**
This implements the master root function as defined in the 
article: https://doi.org/10.1103/PhysRevD.103.023018
**/
CCTK_REAL froot::operator()(const CCTK_REAL mu) 
{
  cache& c{last};
  
  c.lmu               = mu;
  c.x                 = x_from_mu(mu);
  const CCTK_REAL rfsqr  = rfsqr_from_mu_x(mu, c.x);
  const CCTK_REAL qf     = qf_from_mu_x(mu, c.x);
  c.vsqr              = rfsqr * mu*mu;
  
  
  if (c.vsqr >= vsqrinf) {
    c.vsqr = vsqrinf;
    c.w    = winf;
  } else {
    c.w    = 1 / sqrt(1 - c.vsqr);
  }

  c.rho_raw     = d / c.w;
  c.rho         = rho_range.limit_to(c.rho_raw);

  c.eps_raw     = get_eps_raw(mu, qf, rfsqr, c.w);
  c.eps         = eos.range_eps(c.rho, c.ye).limit_to(c.eps_raw); 

  c.press       = eos.at_rho_eps_ye(c.rho, c.eps, c.ye).press();
  ++c.calls;


  const CCTK_REAL a        = c.press / (c.rho * (1. + c.eps));

  const CCTK_REAL h        = (1 + c.eps) * (1 + a);
  
  const CCTK_REAL hbw_raw  = (1 + a) * (1 + qf - mu * rfsqr);
  const CCTK_REAL hbw      = max(hbw_raw, h / c.w);     
  const CCTK_REAL newmu    = 1 / (hbw + rfsqr * mu);
  
  return mu - newmu; 
}


bool froot::stopif(CCTK_REAL mu, CCTK_REAL dmu, CCTK_REAL acc) const
{
  return fabs(dmu) * last.w * last.w < mu * acc;
}


auto froot::initial_bracket(report& errs) const -> interval<CCTK_REAL>
{
  CCTK_REAL mu_max {1.0 / h0};

  if (rsqr >= h0*h0) { 
    const int ndigits2{ 36 };
    const CCTK_REAL margin{ pow(2., 3-ndigits2) };

    f_upper g(h0, rsqr, rbsqr, bsqr);

    ROOTSTAT status;
    mu_max = findroot_using_deriv(g, status, ndigits2, ndigits2+4);

    if (status != ROOTSTAT::SUCCESS) {
      if (status == ROOTSTAT::NOT_CONVERGED) {
        errs.set_prep_root_conv();
      }
      else if (status == ROOTSTAT::NOT_BRACKETED) {
        errs.set_prep_root_bracket();
      }
      return {0, 1.0/h0};
    }
    
    mu_max *= 1. + margin;
    
    assert(g(mu_max).first > 0);
  }
  
  return {0., mu_max};
}


void con2prim_mhd::set_to_nan(prim_vars_mhd& pv, 
                              cons_vars_mhd& cv)
{
  pv.set_to_nan();
  cv.set_to_nan();
}


void con2prim_mhd::operator()(prim_vars_mhd& pv, cons_vars_mhd& cv, 
                               const sm_metric3& g, report& errs) const
{
  errs.iters        = 0;
  errs.adjust_cons  = false;
  errs.set_atmo     = false;
  errs.status       = report::SUCCESS;

  if ((!isfinite(g.vol_elem)) || (g.vol_elem <= 0)) {
    errs.set_invalid_detg(g.vol_elem);
    set_to_nan(pv, cv);
    return;
  } 

  pv.B      = cv.bcons / g.vol_elem;

  const CCTK_REAL d = cv.dens / g.vol_elem;

  if (d <= atmo.rho_cut) {
    errs.set_atmo_set();
    atmo.set(pv, cv, g);
    return;
  }

  const sm_vec3u bu   = cv.bcons / (g.vol_elem * sqrt(d));
  const sm_vec3l rl   = cv.scon / cv.dens;

  const sm_vec3u ru   = g.raise(rl);
  const CCTK_REAL rsqr   = ru * rl;
  const CCTK_REAL rb     = rl * bu;
  const CCTK_REAL rbsqr  = rb * rb;
  const CCTK_REAL bsqr   = g.contract(bu, bu);
  const CCTK_REAL q      = cv.tau / cv.dens;
  const CCTK_REAL ye0    = cv.tracer_ye / cv.dens;

  if ((!isfinite(d)) || (!isfinite(rsqr))  || (!isfinite(ye0)) ||
      (!isfinite(q)) || (!isfinite(rbsqr)) || (!isfinite(bsqr))) 
  {
    errs.set_nans_in_cons(d, q, rsqr, rbsqr, bsqr, ye0);
    set_to_nan(pv, cv);
    return;
  }     

  if (bsqr < 0) 
  {
    errs.set_neg_bsqr(bsqr);
    set_to_nan(pv, cv);
    return;
  }

  if (bsqr > bsqr_lim) 
  {
    errs.set_b_limit(bsqr);
    set_to_nan(pv, cv);
    return; 
  }
  
  const CCTK_REAL ye = eos.range_ye().limit_to(ye0);


  froot::cache sol{};
  froot f{eos, ye, d, q, rsqr, rbsqr, bsqr, sol}; 

  auto bracket = f.initial_bracket(errs);
  
  if (errs.failed()) {
    set_to_nan(pv, cv);
    return;
  }
  
  rarecase nc(bracket, eos.range_rho(), f);

  if (nc.rho_too_big) 
  {
    errs.set_range_rho(d, d);
    set_to_nan(pv, cv);
    return;
  }

  if (nc.rho_too_small) 
  {
    errs.set_atmo_set();
    atmo.set(pv, cv, g);
    return;
  }

  
  ROOTSTAT status;  
  bracket = findroot_no_deriv(f, nc.bracket, acc, max_iter, status);
  
  errs.iters = sol.calls;
  if (status != ROOTSTAT::SUCCESS) {
    if (status == ROOTSTAT::NOT_CONVERGED) {
      errs.set_root_conv();
    }
    else if (status == ROOTSTAT::NOT_BRACKETED) {
      if (nc.rho_big) { 
        errs.set_range_rho(d, d);
        set_to_nan(pv, cv);
        return;
      }
      if (nc.rho_small) {
        errs.set_atmo_set();
        atmo.set(pv, cv, g);
        return;
      }
      errs.set_root_bracket();
    }
    set_to_nan(pv, cv);
    return;
  }
  assert(bracket.contains(sol.lmu));

  
  if (sol.rho < atmo.rho_cut) {
    errs.set_atmo_set();
    atmo.set(pv, cv, g);
    return;
  }

  auto rgeps = eos.range_eps(sol.rho, sol.ye);
  if (sol.eps_raw > rgeps) {
    errs.adjust_cons = true;
    if (sol.rho >= rho_strict) {
      errs.set_range_eps(sol.eps_raw);
      set_to_nan(pv, cv);
      return;
    }
  }
  else if ( sol.eps_raw < rgeps ) {
    errs.adjust_cons = true;    
  }
    
  
  if (! eos.range_ye().contains(ye0) ) {
    errs.adjust_cons = true;
    if ((!ye_lenient) && (sol.rho >= rho_strict)) {
      errs.set_range_ye(ye0);
      set_to_nan(pv, cv);
      return;
    }
  }

  pv.rho    = sol.rho;
  pv.eps    = sol.eps;
  pv.ye     = sol.ye;
  pv.press  = sol.press;
  pv.vel    = sol.lmu * sol.x  * (ru + (rb * sol.lmu) * bu);
  pv.w_lor  = sol.w;

  CCTK_REAL sol_v = sqrt(sol.vsqr);
  if (sol_v > v_lim) {
    pv.rho        = d / w_lim;
    if (pv.rho >= rho_strict) {
      errs.set_speed_limit(sol_v);
      set_to_nan(pv, cv);
      return;
    }
    pv.vel       *= v_lim / sol_v;
    pv.w_lor      = w_lim;

    pv.eps = eos.range_eps(pv.rho, pv.ye).limit_to(pv.eps); 

    pv.press       = eos.at_rho_eps_ye(pv.rho, pv.eps, pv.ye).press();
    
    errs.adjust_cons = true;   
  }

  sm_vec3l El = g.cross_product(pv.B, pv.vel);
  pv.E = g.raise(El);

  if (errs.adjust_cons) {
    cv.from_prim(pv, g);
  }
}





f_rare::f_rare(CCTK_REAL wtarg_, const froot& f_)
: v2targ(1.0 - 1.0/(wtarg_*wtarg_)), f(f_) {}


/**
This implements a root function for finding mu from W_hat
**/
auto f_rare::operator()(const CCTK_REAL mu) const 
-> std::pair<CCTK_REAL, CCTK_REAL>
{
  CCTK_REAL x     = f.x_from_mu(mu);
  CCTK_REAL xsqr  = x*x;
  CCTK_REAL rfsqr = f.rfsqr_from_mu_x(mu,x);
  CCTK_REAL vsqr  = mu*mu * rfsqr;  
  
  CCTK_REAL y     = vsqr - v2targ;
  CCTK_REAL dy = 2.0*mu*x*(xsqr * f.rsqr + mu * (xsqr+x+1.0) * f.rbsqr);
  
  return {y, dy};
}


/**
This checks for the corner case where the density might cross the 
allowed bounds while finding the root of the master function. In
that case, a tighter initial root finding interval is constructed
to guarantee uniqueness.
**/
rarecase::rarecase(const interval<CCTK_REAL> ibracket,
                   const interval<CCTK_REAL> rgrho, const froot& f)
{
  
  CCTK_REAL muc0 = ibracket.min();
  CCTK_REAL muc1 = ibracket.max();
  const int ndigits = 30; 
  //We just assume 30 binary digits for finding adjusted interval should 
  //be more than enough. Speed is not an issue for this rare case.
          
  
  if (f.d > rgrho.max()) {
    CCTK_REAL wc = f.d / rgrho.max(); 
    if (wc  > f.winf) {
      rho_too_big = true; 
    }
    else {
      f_rare g(wc, f);
      
      if (g(muc1).first <= 0) {
        rho_too_big = true; 
      }
      else {
        if (g(muc0).first < 0) {
          ROOTSTAT status;
          CCTK_REAL mucc = findroot_using_deriv(g, ibracket, 
                                             status, ndigits, ndigits+2);
          assert(status == ROOTSTAT::SUCCESS); 
          muc0 = max(muc0, mucc);
          rho_big = true;
        }
      }
    }
  }

  if (f.d <  f.winf * rgrho.min()) {
    CCTK_REAL wc = f.d / rgrho.min(); 
    if (wc < 1) {
      rho_too_small = true; 
    }
    else {
      f_rare g(wc, f);
      if (g(muc0).first >= 0) {
        rho_too_small = true; 
      }
      else {
        if (g(muc1).first > 0) {
          ROOTSTAT status;
          CCTK_REAL mucc = findroot_using_deriv(g, ibracket, 
                                             status, ndigits, ndigits+2);
          assert(status == ROOTSTAT::SUCCESS); 
          muc1 = min(muc1, mucc);
          rho_small = true;
        }
      }
    }
  }

  bracket = interval<CCTK_REAL>{muc0, muc1};
}

#endif
