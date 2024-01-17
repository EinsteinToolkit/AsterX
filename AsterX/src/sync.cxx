#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace AsterX {

////////////////////////////////////////////////////////////////////////////////

extern "C" void AsterX_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Sync;

  // do nothing
}

} // namespace AsterX
