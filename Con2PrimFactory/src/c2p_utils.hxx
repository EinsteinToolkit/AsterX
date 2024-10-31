#ifndef C2PUTILS_HXX
#define C2PUTILS_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <mat.hxx>
#include <vec.hxx>
#include <sum.hxx>
#include <simd.hxx>

#include <algorithm>
#include <array>
#include <cmath>

#include "aster_utils.hxx"

namespace Con2PrimFactory {

enum class ROOTSTAT {
  SUCCESS,
  NOT_CONVERGED,
  NOT_BRACKETED
};

} // namespace Con2PrimFactory

#endif // #ifndef C2PUTILS_HXX
