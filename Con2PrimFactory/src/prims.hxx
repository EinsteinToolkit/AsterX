/*! \file prims.hxx
\brief Class definitions representing primitive mhd variables.
*/

#ifndef PRIMS_HXX
#define PRIMS_HXX

#include "c2p_utils.hxx"

namespace Con2PrimFactory {

using namespace AsterUtils;

struct prim_vars {
  CCTK_REAL rho;
  CCTK_REAL eps;
  CCTK_REAL Ye;
  CCTK_REAL press;
  vec<CCTK_REAL, 3> vel;
  CCTK_REAL w_lor;
  vec<CCTK_REAL, 3> Bvec;
  vec<CCTK_REAL, 3> E;

  /// Default constructor. Leaves all members uninitialized.
  prim_vars() = default;

  /// Construct from single variables.
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline prim_vars(
      CCTK_REAL rho_, CCTK_REAL eps_, CCTK_REAL ye_, CCTK_REAL press_,
      vec<CCTK_REAL, 3> vel_, CCTK_REAL w_lor_, vec<CCTK_REAL, 3> Bvec_)
      : rho(rho_), eps(eps_), Ye(ye_), press(press_), vel(vel_), w_lor(w_lor_),
        Bvec(Bvec_){};

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  scatter(CCTK_REAL &rho_, CCTK_REAL &eps_, CCTK_REAL &ye_, CCTK_REAL &press_,
          CCTK_REAL &velx_, CCTK_REAL &vely_, CCTK_REAL &velz_,
          CCTK_REAL &w_lor_, CCTK_REAL &Bvec_x_, CCTK_REAL &Bvec_y_,
          CCTK_REAL &Bvec_z_, CCTK_REAL &E_x_, CCTK_REAL &E_y_,
          CCTK_REAL &E_z_) const {
    rho_ = rho;
    eps_ = eps;
    ye_ = Ye;
    press_ = press;
    velx_ = vel(0);
    vely_ = vel(1);
    velz_ = vel(2);
    w_lor_ = w_lor;
    Bvec_x_ = Bvec(0);
    Bvec_y_ = Bvec(1);
    Bvec_z_ = Bvec(2);
    E_x_ = E(0);
    E_y_ = E(1);
    E_z_ = E(2);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_to_nan() {
    rho = eps = Ye = press = vel(0) = vel(1) = vel(2) = w_lor = Bvec(0) =
        Bvec(1) = Bvec(2) = E(0) = E(1) = E(2) =
            std::numeric_limits<CCTK_REAL>::quiet_NaN();
  }
};

} // namespace Con2PrimFactory
#endif
