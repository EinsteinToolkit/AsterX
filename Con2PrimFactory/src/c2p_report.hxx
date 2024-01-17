// Error report adapted from the RePrimAnd library

#ifndef C2P_REPORT_HXX
#define C2P_REPORT_HXX

namespace Con2PrimFactory {

/// Struct to represent the outcome of a con2prim call.
/**
This should be checked after con2prim. The various members indicate
success or failure, if the conserved variables were changed, or
artificial atmosphere applied.
*/
struct c2p_report {
public:
  /// Failure flags
  enum err_code {
    SUCCESS,                /*!<Recovery succeeded (possibly after
                                allowed corrections) */
    INVALID_DETG,           ///< 3-metric determinant not positive or NAN/INF
    NEG_BSQR,               ///< Negative B square, metric not positive definite
    NANS_IN_CONS,           ///< One or more evolved variables is NAN/INF
    RANGE_RHO,              ///< Mass density outside EOS range
    RANGE_EPS,              ///< Fluid internal energy outside EOS range
    SPEED_LIMIT,            ///< Speed limit exceeded
    RANGE_YE,               ///< Electron fraction outside EOS range
    B_LIMIT,                ///< Magnetization above limit
    ROOT_FAIL_CONV,         /*!<Root finding did not converge (this
                                should never happen) */
    ROOT_FAIL_BRACKET,      /*!<Bracketing of root failed (this
                                should never happen) */
    PREP_ROOT_FAIL_CONV,    /*!<Auxiliary root finding did not
                                converge (this should never happen) */
    PREP_ROOT_FAIL_BRACKET, /*!<Bracketing of auxiliary root
                                failed (this should never happen) */
    ERR_CODE_NOT_SET
  };

  /// Default constructor, resulting object invalid.
  c2p_report() = default;

  /**
  @return String with debug information.

  This assembles a human readable debug message from the flags
  and numerical values stored internally (for performance reasons, no
  string is ever created before explicitly requested via this function).
  **/
  //  void debug_message() const;

  /**
  @return If the input was invalid according to the error policy.
  **/
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  failed() const {
    return status != SUCCESS;
  }

  /// SUCCESS or reason for failure.
  err_code status{ERR_CODE_NOT_SET};

  /// Whether the conserved variables were adjusted.
  bool adjust_cons{false};

  /// Whether the artificial atmosphere was enforced.
  bool set_atmo{false};

  /// Number of calls to the EOS needed for the root finding.
  CCTK_INT iters;

private:
  /// Conserved density \f$ D \f$.
  CCTK_REAL dens;

  /// Square of conserved momentum \f$ S^2 = { S_i S^i} \f$.
  CCTK_REAL Ssq;

  /// Square of B \f$ B^2 = B^i B_i \f$.
  CCTK_REAL Bsq;

  /// Product \f$ (B^i S_i)^2 \f$.
  CCTK_REAL BiSi;

  // /// Specific total energy density \f$ q = \frac{\tau}{D}  \f$.
  // CCTK_REAL qtot;

  /// Electron fraction \f$ Y_e \f$.
  CCTK_REAL Ye;

  /// Rest mass density \f$ \rho \f$.
  CCTK_REAL rho;

  /// Specific internal energy \f$ \epsilon \f$.
  CCTK_REAL eps;

  /// Velocity.
  vec<CCTK_REAL, 3> vel;

  /// 3-Metric determinant.
  CCTK_REAL detg;

  // CXX implementation

public:
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_atmo_set() {
    status = SUCCESS;
    set_atmo = true;
    adjust_cons = true;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_invalid_detg(CCTK_REAL detg_) {
    status = INVALID_DETG;
    adjust_cons = true;
    detg = detg_;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_neg_Bsq(CCTK_REAL Bsq_) {
    status = NEG_BSQR;
    adjust_cons = true;
    Bsq = Bsq_;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_nans_in_cons(CCTK_REAL dens_, CCTK_REAL Ssq_, CCTK_REAL Bsq_,
                   CCTK_REAL BiSi_, CCTK_REAL Ye_) {
    status = NANS_IN_CONS;
    set_atmo = false;
    adjust_cons = true;
    dens = dens_;
    Ssq = Ssq_;
    Bsq = Bsq_;
    BiSi = BiSi_;
    Ye = Ye_;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_range_rho(CCTK_REAL dens_, CCTK_REAL rho_) {
    status = RANGE_RHO;
    set_atmo = false;
    adjust_cons = true;
    dens = dens_;
    rho = rho_;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_range_eps(CCTK_REAL eps_) {
    status = RANGE_EPS;
    set_atmo = false;
    adjust_cons = false; // we do not adjust cons in this case!
    eps = eps_;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_speed_limit(vec<CCTK_REAL, 3> vel_) {
    status = SPEED_LIMIT;
    set_atmo = false;
    adjust_cons = true;
    vel = vel_;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_B_limit(CCTK_REAL Bsq_) {
    status = B_LIMIT;
    set_atmo = false;
    adjust_cons = true;
    Bsq = Bsq_;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_range_Ye(CCTK_REAL Ye_) {
    status = RANGE_YE;
    set_atmo = false;
    adjust_cons = true;
    Ye = Ye_;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_root_conv() {
    status = ROOT_FAIL_CONV;
    set_atmo = false;
    adjust_cons = true;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_root_bracket() {
    status = ROOT_FAIL_BRACKET;
    set_atmo = false;
    adjust_cons = true;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_prep_root_conv() {
    status = PREP_ROOT_FAIL_CONV;
    set_atmo = false;
    adjust_cons = true;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_prep_root_bracket() {
    status = PREP_ROOT_FAIL_BRACKET;
    set_atmo = false;
    adjust_cons = true;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  debug_message() const {
    switch (status) {
    case SUCCESS:
      printf("Con2Prim succeeded. \n");
      if (set_atmo) {
        printf("Artificial atmosphere has been enforced. \n");
      }
      if (adjust_cons) {
        printf("Conserved variables have been changed. \n");
      }
      break;
    case INVALID_DETG:
      printf("Invalid metric determinant: detg = %16.8e \n", detg);
      break;
    case NEG_BSQR:
      printf("3-metric not positive, negative B^2: Bsq = %16.8e \n", Bsq);
      break;
    case NANS_IN_CONS:
      printf("NAN in conserved variables,"
             "dens = %16.8e, Ssq = %16.8e, Bsq = %16.8e, BiSi = %16.8e,"
             "Ye = %16.8e \n",
             dens, Ssq, Bsq, BiSi, Ye);
      break;
    case RANGE_RHO:
      printf("Density out of range, dens = %16.8e, rho = %16.8e \n", dens, rho);
      break;
    case RANGE_EPS:
      printf("Specific energy out of range, eps = %16.8e \n", eps);
      break;
    case SPEED_LIMIT:
      printf("Speed limit exceeded, vx, vy, vz = %16.8e, %16.8e, %16.8e \n",
             vel(0), vel(1), vel(2));
      break;
    case B_LIMIT:
      printf("Limit for magnetic field exceeded, B = %16.8e \n", sqrt(Bsq));
      break;
    case RANGE_YE:
      printf("Electron fraction out of range, Y_e = %16.8e \n", Ye);
      break;
    case ROOT_FAIL_CONV:
      printf("Root finding failed (not converged after %i steps) \n", iters);
      break;
    case ROOT_FAIL_BRACKET:
      printf("Root finding failed (faulty bracketing) \n");
      break;
    case PREP_ROOT_FAIL_CONV:
      printf("Preparatory root finding failed (not converged) \n");
      break;
    case PREP_ROOT_FAIL_BRACKET:
      printf("Preparatory root finding failed (faulty bracketing) \n");
      break;
    default:
      assert(false);
      printf("Invalid error type. Should never happen. Code is messed up. \n");
      break;
    }
  }

}; // struct

} // namespace Con2PrimFactory
#endif
