#ifndef C2P_1DREPRIMAND_ROOTFINDER_HXX
#define C2P_1DREPRIMAND_ROOTFINDER_HXX

#include <boost/math/tools/roots.hpp> //TODO: fix dependency

#include "c2p.hxx"
#include "roots.hxx"
#include "c2p_1DRePrimAnd_intervals.hxx"

namespace Con2PrimFactory {

template<class F, class T = typename F::value_t> 
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline auto findroot_using_deriv(F& f, 
      interval<T> bracket, ROOTSTAT& errs, 
      int digits, unsigned int max_calls=20) -> T
{
  if (max_calls < 4) {
    printf("c2p_1DRePrimAnd_rootfinder: Root finding call limit set too low for "
                           "meaningful results");
  }
  
  T f_left  { f(bracket.min()).first }; 
  T f_right { f(bracket.max()).first }; 
    
  max_calls -= 2;
    
  if (f_left * f_right >= 0) {
    if (f_right == 0) {
      errs = ROOTSTAT::SUCCESS;
      return bracket.max();
    }
    if (f_left == 0) {
      errs = ROOTSTAT::SUCCESS;
      return bracket.min();
    }    
    errs = ROOTSTAT::NOT_BRACKETED;
    return std::numeric_limits<T>::quiet_NaN();    
  }


  T x_guess { (bracket.min() * f_right - bracket.max() * f_left) 
                / (f_right - f_left) };

  boost::uintmax_t iters {max_calls};
  T sol { boost::math::tools::newton_raphson_iterate(f, x_guess, 
                    bracket.min(), bracket.max(), digits, iters) };

  errs = (iters == max_calls) 
         ? ROOTSTAT::NOT_CONVERGED : ROOTSTAT::SUCCESS;

  return sol;  
}

template<class F, class T = typename F::value_t> 
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline auto findroot_using_deriv(F& f, ROOTSTAT& errs, 
      int digits, unsigned int max_calls=20) -> T
{
  return findroot_using_deriv(f, f.initial_bracket(), 
                              errs, digits, max_calls);
}


template<class F, class T = typename  F::value_t> 
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline auto findroot_no_deriv(F& f, interval<T> bracket, T tol, 
       unsigned int max_calls, ROOTSTAT& errs) -> interval<T>
{
  if (max_calls < 10) {
    printf("c2p_1DRePrimAnd_rootfinder: Root finding call limit set too low for "
                      "meaningful results");
  }
  
  T f_left  { f(bracket.min()) }; 
  T f_right { f(bracket.max()) }; 
    
  if (f_left * f_right >= 0) {
    if (f_right==0) {
      errs = ROOTSTAT::SUCCESS;
      return {bracket.max(),bracket.max()};
    }
    if (f_left==0) {
      errs = ROOTSTAT::SUCCESS;
      return {bracket.min(),bracket.min()};
    }
    errs = ROOTSTAT::NOT_BRACKETED;
    return {std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::max()};    
  }
  
  auto stopif = [&] (T l, T r) {
    return f.stopif(l, r-l, tol);
  };
  
  max_calls    -= 2;
  boost::uintmax_t iters{ max_calls };
  
  auto res = boost::math::tools::toms748_solve(
    f, bracket.min(), bracket.max(), f_left, f_right,
    stopif, iters
  );
  
  errs = (iters == max_calls) 
         ? ROOTSTAT::NOT_CONVERGED : ROOTSTAT::SUCCESS;
  
  return {res.first, res.second};
}

template<class F, class T = typename  F::value_t> 
auto findroot_no_deriv(F& f, T tol, 
       unsigned int max_calls, ROOTSTAT& errs) -> interval<T>
{
  return findroot_no_deriv(f, f.initial_bracket(), 
                           tol, max_calls, errs);
}


}

#endif
