#include <loop.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C" void AsterMasks_Initialize(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTSX_AsterMasks_Initialize;
  DECLARE_CCTK_PARAMETERS;
 
  grid.loop_all_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
  
          aster_mask_vc(p.I) = 1.0;

  });
}


