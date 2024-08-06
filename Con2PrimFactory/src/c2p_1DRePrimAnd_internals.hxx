/*! \file c2p_1DRePrimAnd_internals.hxx
\brief Internal stuff for primitive recovery algorithm.
*/

#ifndef C2P_1DREPRIMAND_INTERNALS_HXX
#define C2P_1DREPRIMAND_INTERNALS_HXX

#include "c2p.hxx"

namespace Con2PrimFactory {

class rarecase;
class f_rare;

/// Function object representing the root function.
/** This contains all the fixed parameters defining the function.
    It also remembers intermediate results from the last evaluation,
    which we do not want to recompute for performance reasons.
    In particular, it remembers if specific internal energy and density
    were in the validity region of the EOS, or if they needed to be limited.
**/
class froot {
  using range   = eos::range;
  using report  = c2p_report;

  //const eos_thermal eos;    ///< The EOS.
  const CCTK_REAL h0;         ///< Lower bound for enthalpy, \f$ h_0 \f$
  const range rho_range;    ///< Valid density interval of the EOS. 
  const CCTK_REAL d;          ///< \f$ d = \frac{D}{\sqrt{\det(g_{ij})}} \f$
  const CCTK_REAL qtot;       ///< \f$ q = \frac{\tau}{D}  \f$
  const CCTK_REAL rsqr;       ///< \f$ r^2 = \frac{ S_i S^i}{D^2} \f$
  const CCTK_REAL rbsqr;      ///< \f$ (r^l b_l)^2 \f$
  const CCTK_REAL bsqr;       ///< \f$ b^2 = \frac{B^2}{D} \f$
  const CCTK_REAL brosqr;     ///< \f$ b^2 r^2_\perp = b^2 r^2 - (r^lb_l)^2 \f$
  CCTK_REAL winf;             ///< Upper bound for Lorentz factor
  CCTK_REAL vsqrinf;          ///< Upper bound for squared velocity 
  
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL x_from_mu(CCTK_REAL mu) const;
  
  ///Computes fluid momentum from total one, for a given \f$ \mu, x \f$.
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL rfsqr_from_mu_x(CCTK_REAL mu, CCTK_REAL x) const;
  
  ///Computes fluid energy from total one, for a given \f$ \mu, x \f$.
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL qf_from_mu_x(CCTK_REAL mu, CCTK_REAL x) const;
  
  ///Computes specific energy from conserved fluid variables and velocity. 
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline static CCTK_REAL get_eps_raw(CCTK_REAL mu,  CCTK_REAL qf, 
                            CCTK_REAL rfsqr, CCTK_REAL w);  
  
  friend class rarecase;
  friend class f_rare;
     
  public:
  
  using value_t = CCTK_REAL; 


  
  ///Store intermediate results obtained when evaluating function.
  struct cache {
    CCTK_REAL ye;       ///< Electron fraction in valid range.
    CCTK_REAL lmu;      ///< \f$ \mu = \frac{1}{h W} \f$.  
    CCTK_REAL x;         ///< \f$ \frac{1}{1 + \mu b^2} \f$.
    CCTK_REAL rho;      ///< Rest mass density \f$ \rho \f$.
    CCTK_REAL rho_raw;  ///< Rest mass density not limited to EOS range.
    CCTK_REAL eps;      ///< Specific internal energy, limited to EOS range.
    CCTK_REAL eps_raw;  ///< Specific internal energy, not limited.  
    CCTK_REAL press;    ///< Pressure \f$ P \f$.
    CCTK_REAL vsqr;     ///< Squared 3-velocity \f$ v^2 \f$.  
    CCTK_REAL w;         ///< Lorentz factor \f$ W \f$.
    unsigned int calls;
  };
  
  /// Constructor
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline froot(
    // const eos_thermal& eos_,       ///< The EOS
    CCTK_REAL valid_ye,               ///< Electron fraction
    CCTK_REAL d_,                     ///< \f$ d = \frac{D}{\sqrt{\det(g_{ij})}} \f$
    CCTK_REAL qtot_,                  ///< \f$ q = \frac{\tau}{D} + EM\f$
    CCTK_REAL rsqr_,                  ///< \f$ r^2 = \frac{ S_i S^i}{D^2} \f$
    CCTK_REAL rbsqr_,                 ///< \f$ (r^l b_l)^2 \f$
    CCTK_REAL bsqr_,                  ///< \f$ b^2 = \frac{B^2}{D} \f$
    cache& last_                   ///< cache for intermediate results
  );

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline froot(const froot&)            = default;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline froot(froot&&)                 = default;


  /// The root function
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL operator()(CCTK_REAL mu);

  /// The convergence criterion for root finding.
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool stopif(CCTK_REAL mu, CCTK_REAL dmu, CCTK_REAL acc) const;

  /// Initial guess for root finding
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline auto initial_bracket(report& errs) const -> interval<CCTK_REAL>;
  
  private:
  
  cache& last;


};

///Class representing the auxiliary root function 
class f_upper {
  public:
  using value_t = CCTK_REAL;
  
  /// Constructor
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline f_upper(
    CCTK_REAL h0_,         ///< Lower bound for enthalpy
    CCTK_REAL rsqr_,       ///< \f$ r^2 = \frac{ S_i S^i}{D^2} \f$
    CCTK_REAL rbsqr_,      ///< \f$ (r^l b_l)^2 \f$
    CCTK_REAL bsqr_        ///< \f$ b^2 = \frac{B^2}{D} \f$
  );

  
  /// The function and first derivative
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline auto operator()(CCTK_REAL mu) const -> std::pair<CCTK_REAL,CCTK_REAL>;


  /// Initial bracket for root finding
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline auto initial_bracket() const -> interval<CCTK_REAL>;
  
  private:
  
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline auto x_from_mu(CCTK_REAL mu) const -> CCTK_REAL;
  
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline auto rfsqr_from_mu_x(CCTK_REAL mu, CCTK_REAL x) const -> CCTK_REAL;
  
  /// Compute new \f$ h_0 W \f$ from initial guess.
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline auto new_h0w_from_mu_x(CCTK_REAL mu, CCTK_REAL x) const -> CCTK_REAL;
  
  const CCTK_REAL h0;      ///< Lower bound for enthalpy, \f$ h_0 \f$
  const CCTK_REAL h0sqr;   ///< \f$ h_0^2 \f$
  const CCTK_REAL rsqr;    ///< Fixed parameter \f$ r^2 = \frac{ S_i S^i}{D^2} \f$
  const CCTK_REAL rbsqr;   ///< Fixed parameter \f$ (r^l b_l)^2 \f$
  const CCTK_REAL bsqr;    ///< Fixed parameter \f$ b^2 = \frac{B^2}{D} \f$
};




///Root function used by rarecase class   
class f_rare {
  const CCTK_REAL v2targ;  ///< Target squared velocity   
  const froot& f;

  public:
  using value_t = CCTK_REAL;
  
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline f_rare(
    CCTK_REAL wtarg_,      ///< Target Lorentz factor
    const froot& f_
  );

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline auto operator()(CCTK_REAL mu) const -> std::pair<CCTK_REAL,CCTK_REAL>;  
};

///Class for handling rare corner case 
class rarecase {  
  public:
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline rarecase(
    const interval<CCTK_REAL> bracket,   ///< Initial master root bracket 
    const interval<CCTK_REAL> rgrho,     ///< Allowed density range
    const froot& f                    ///< Master root function
  );

  /// Root bracket on which solution is unique
  interval<CCTK_REAL> bracket;  
  
  bool rho_too_big { false };    ///< Density definitely too large
  bool rho_big { false };        ///< Possibly too large
  bool rho_too_small { false };  ///< Density definitely too small
  bool rho_small { false };      ///< Possibly too small
};

}

#endif
