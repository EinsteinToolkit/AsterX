#ifndef C2P_1DREPRIMAND_INTERVALS_HXX
#define C2P_1DREPRIMAND_INTERVALS_HXX

#include <cmath>
#include <algorithm>
#include <cassert>
#include <loop_device.hxx>

namespace Con2PrimFactory {

/**Class representing a closed interval

@tparam T the underlying numerical type, e.g. double
**/
template<class T>
class interval {
  T min_{0};  ///< Lower boundary 
  T max_{0};  ///< Upper boundary
  public:
  
  ///Default constructor yields empty range \f$ [0,0] \f$.
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline interval() noexcept = default;
  
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline ~interval() noexcept = default;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline interval(const interval&) noexcept = default;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline interval(interval&&) noexcept = default;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline interval& operator=(const interval&) noexcept = default;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline interval& operator=(interval&&) noexcept = default;
  
  /**\brief Construct from minimum and maximum
  
  @param min Lower boundary of interval
  @param max Upper boundary of interval
  
  \pre Aborts unless max >= min 
  \pre Aborts unless max and min are both finite values
  **/
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline interval(T min, T max) : min_(min), max_(max) 
  {
    assert(min_ <= max_);
    assert(std::isfinite(min_));
    assert(std::isfinite(max_));
  }
  
  /**
  @param x Value to test
  @return If value is contained in closed interval \f$[a,b]\f$. 
  
  \note Returns false if x is NAN or INF.
  **/
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool contains(const T& x) const 
  {
    return (x >= min_) && (x <= max_);
  }

  /**
  @param r Interval to test
  @return If Interval is contained in closed interval \f$[a,b]\f$.   
  **/

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool contains(const interval<T>& r) const 
  {
    return (contains(r.min()) && contains(r.max()));
  }
  
  /**\brief Limit value to interval
  
  @param x Value to constrain
  @return x if it inside interval, else the closest boundary
  **/
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T limit_to(const T& x) const 
  {
    return std::min(std::max(min_, x), max_);
  }
  
  /**@return Lower boundary **/
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T min() const noexcept {return min_;}
  
  /**@return Upper boundary **/
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T max() const noexcept {return max_;}
  
  /** @return Length of interval **/
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T length() const {return max_ - min_;}
};

///Test if value is above interval
template<class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool operator>(T x, const interval<T>& i) 
{
  return x > i.max();  
}

template<class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool operator>=(T x, const interval<T>& i) 
{
  return x >= i.max();  
}

///Test if value is below interval
template<class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool operator<(T x, const interval<T>& i) 
{
  return x < i.min();
}

template<class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool operator<=(T x, const interval<T>& i) 
{
  return x <= i.min();
}


template<class T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline auto intersect(const interval<T>& a, const interval<T>& b)
-> interval<T>
{
  return {std::max(a.min(), b.min()), std::min(a.max(), b.max())};
}

template<class T, class... C>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline auto intersect(const interval<T>& a, const interval<T>& b, 
               const C&... c)
-> interval<T>
{
  return intersect(intersect(a,b), c...);
}



}// namespace Con2PrimFactory

#endif
