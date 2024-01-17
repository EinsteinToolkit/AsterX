#ifndef EOS_UTILS_HXX
#define EOS_UTILS_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace EOSX {

using namespace std;
using namespace Loop;
using namespace Arith;

/// Class representing a range
struct eos_range {
  CCTK_REAL min; ///< Minimum
  CCTK_REAL max; ///< Maximum
  /// Default constructor: empty range.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_range()
      : min(0), max(0) {}
  /// Construct from minimum and maximum
  CCTK_DEVICE
  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_range(CCTK_REAL min_,
                                                          CCTK_REAL max_)
      : min(min_), max(max_) {}
  /// Check if value is contained in [min,max]. False for NAN or INF.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  contains(CCTK_REAL x) const {
    return (x >= min) && (x <= max);
  }
};

// TODO: add enums for error messages
/// Class representing error conditions in EOS calls.
struct eos_status {
  bool failed; ///< Set to true if parameters are out of range/NAN/INF
  //  std::string  err_msg; ///< Error description in case of failure, else
  //  undefined.
  /// Default constructor: Set to no failure.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_status()
      : failed(false) {}
  /// Set fail flag and error message.
  //  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  //  fail(std::string msg) {
  //    failed = true;
  //    err_msg = msg;
  //  }
};

} // namespace EOSX

#endif // #ifndef UTILS_HXX
