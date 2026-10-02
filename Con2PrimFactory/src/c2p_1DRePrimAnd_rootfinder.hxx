#ifndef C2P_1DREPRIMAND_ROOTFINDER_HXX
#define C2P_1DREPRIMAND_ROOTFINDER_HXX

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

#include "roots.hxx" // Algo::brent
#include "c2p_1DRePrimAnd_intervals.hxx"
#include "c2p_toms748.hxx"
#include "c2p_utils.hxx"

namespace Con2PrimFactory {

namespace RePrimAnd_rootdetail {

template <class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool same_sign(T a,
                                                                        T b) {
  return (a > T(0) && b > T(0)) || (a < T(0) && b < T(0));
}

template <class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
in_open_interval(T x, T a, T b) {
  return (x > a) && (x < b);
}

template <class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T midpoint(T a,
                                                                     T b) {
  return (a + b) * T(0.5);
}

template <class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T secant(T a, T fa,
                                                                   T b, T fb) {
  const T denom = (fb - fa);
  if (denom == T(0))
    return std::numeric_limits<T>::quiet_NaN();
  return (a * fb - b * fa) / denom;
}

} // namespace RePrimAnd_rootdetail

// Derivative-assisted Newton with bracketing, Brent fallback on pathology.
template <class F, class T = typename F::value_t>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
findroot_using_deriv(F &f, interval<T> bracket, ROOTSTAT &errs, int digits,
                     unsigned int max_calls = 256) {
  using namespace RePrimAnd_rootdetail;

  if (max_calls < 4u) {
    errs = ROOTSTAT::NOT_CONVERGED;
    return std::numeric_limits<T>::quiet_NaN();
  }

  T a = bracket.min();
  T b = bracket.max();
  if (b < a) {
    const T tmp = a;
    a = b;
    b = tmp;
  }

  auto eval_pair = [&](T x) -> std::pair<T, T> {
    const auto val = f(x);
    return {T(val.first), T(val.second)};
  };

  auto va = eval_pair(a);
  auto vb = eval_pair(b);
  T fa = va.first;
  T fb = vb.first;

  // Match RePrimAnd endpoint precedence: right endpoint first.
  if (fb == T(0)) {
    errs = ROOTSTAT::SUCCESS;
    return b;
  }
  if (fa == T(0)) {
    errs = ROOTSTAT::SUCCESS;
    return a;
  }
  if (!std::isfinite(fa) || !std::isfinite(fb) || same_sign(fa, fb)) {
    errs = ROOTSTAT::NOT_BRACKETED;
    return std::numeric_limits<T>::quiet_NaN();
  }

  const int minbits =
      (digits > 0) ? digits : (std::numeric_limits<T>::digits - 4);
  const T tol = std::ldexp(T(1), -minbits);
  const unsigned int newton_calls = max_calls - 2u;

  // Same initial guess as RePrimAnd reference.
  T x = secant(a, fa, b, fb);
  if (!std::isfinite(x) || !in_open_interval(x, a, b))
    x = midpoint(a, b);

  int stalled = 0;
  for (unsigned int calls = 0; calls < newton_calls; ++calls) {
    const auto vx = eval_pair(x);
    const T fx = vx.first;
    const T dfx = vx.second;

    if (!std::isfinite(fx) || !std::isfinite(dfx))
      break;
    if (fx == T(0)) {
      errs = ROOTSTAT::SUCCESS;
      return x;
    }

    if (same_sign(fa, fx)) {
      a = x;
      fa = fx;
    } else {
      b = x;
      fb = fx;
    }

    const T width = b - a;
    const T scale = std::max(T(1), std::max(std::abs(a), std::abs(b)));
    if (std::abs(width) <= tol * scale) {
      errs = ROOTSTAT::SUCCESS;
      return x;
    }

    T x_new = std::numeric_limits<T>::quiet_NaN();
    const T dfx_floor = std::sqrt(std::numeric_limits<T>::epsilon());
    if (std::abs(dfx) > dfx_floor)
      x_new = x - fx / dfx;
    if (!std::isfinite(x_new) || !in_open_interval(x_new, a, b))
      x_new = secant(a, fa, b, fb);
    if (!std::isfinite(x_new) || !in_open_interval(x_new, a, b))
      x_new = midpoint(a, b);

    const T step_scale = std::max(T(1), std::abs(x));
    if (std::abs(x_new - x) <= std::numeric_limits<T>::epsilon() * step_scale)
      ++stalled;
    else
      stalled = 0;
    x = x_new;

    if (stalled >= 2)
      break;
  }

  // Fallback to robust bracketing in pathological derivative/stall cases.
  auto fn = [&](T x0) -> T { return eval_pair(x0).first; };
  if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(fa) ||
      !std::isfinite(fb) || same_sign(fa, fb)) {
    errs = ROOTSTAT::NOT_BRACKETED;
    return std::numeric_limits<T>::quiet_NaN();
  }

  int iters = 0;
  bool failed = false;
  auto ab = Algo::brent(fn, a, b, minbits, int(max_calls), iters, failed);
  errs = failed ? ROOTSTAT::NOT_CONVERGED : ROOTSTAT::SUCCESS;
  if (errs != ROOTSTAT::SUCCESS)
    return std::numeric_limits<T>::quiet_NaN();

  const T a_root = ab.first, b_root = ab.second;
  const T fa2 = fn(a_root), fb2 = fn(b_root);
  if (!std::isfinite(fa2) || !std::isfinite(fb2)) {
    errs = ROOTSTAT::NOT_CONVERGED;
    return std::numeric_limits<T>::quiet_NaN();
  }

  if (fb2 == T(0) || std::abs(fb2) < std::abs(fa2))
    return b_root;
  if (std::abs(fa2) < std::abs(fb2))
    return a_root;
  return (a_root + b_root) * T(0.5);
}

// Derivative-free TOMS748 adapter, Brent fallback on pathological stall.
template <class F, class T = typename F::value_t>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline interval<T>
findroot_no_deriv(F &f, interval<T> bracket, T acc, unsigned int max_calls,
                  ROOTSTAT &errs) {
  if (max_calls < 10u) {
    errs = ROOTSTAT::NOT_CONVERGED;
    T a = bracket.min();
    T b = bracket.max();
    if (b < a) {
      const T tmp = a;
      a = b;
      b = tmp;
    }
    return {a, b};
  }

  interval<T> res = toms748_solve(f, bracket, acc, max_calls, errs);
  if (errs == ROOTSTAT::SUCCESS || errs == ROOTSTAT::NOT_BRACKETED)
    return res;

  T a = res.min();
  T b = res.max();
  if (b < a) {
    const T tmp = a;
    a = b;
    b = tmp;
  }
  if (!std::isfinite(a) || !std::isfinite(b))
    return res;

  auto fn = [&](T x) -> T { return f(x); };
  const T fa = fn(a), fb = fn(b);
  if (fb == T(0))
    return {b, b};
  if (fa == T(0))
    return {a, a};
  if (!std::isfinite(fa) || !std::isfinite(fb) || fa * fb > T(0)) {
    errs = ROOTSTAT::NOT_BRACKETED;
    return {std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max()};
  }

  const int minbits = std::numeric_limits<T>::digits - 4;
  int iters = 0;
  bool failed = false;
  auto ab = Algo::brent(fn, a, b, minbits, int(max_calls), iters, failed);
  errs = failed ? ROOTSTAT::NOT_CONVERGED : ROOTSTAT::SUCCESS;

  return {ab.first, ab.second};
}

template <class F, class T = typename F::value_t>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
findroot_using_deriv(F &f, ROOTSTAT &errs, int digits,
                     unsigned int max_calls = 256) {
  return findroot_using_deriv(f, f.initial_bracket(), errs, digits, max_calls);
}

template <class F, class T = typename F::value_t>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline interval<T>
findroot_no_deriv(F &f, T tol, unsigned int max_calls, ROOTSTAT &errs) {
  return findroot_no_deriv(f, f.initial_bracket(), tol, max_calls, errs);
}

} // namespace Con2PrimFactory
#endif
