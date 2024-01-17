#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>


namespace KadathImporterX {
extern "C"
void KadathImporterX_check_parameters (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("Checking ID parameters");

  if (not CCTK_EQUALS (initial_data,    "Kadath") or
      not CCTK_EQUALS (initial_lapse,   "Kadath") or
      not CCTK_EQUALS (initial_shift,   "Kadath") or
      not (CCTK_EQUALS (initial_dtlapse, "Kadath") or
           CCTK_EQUALS (initial_dtlapse, "none") or
           CCTK_EQUALS (initial_dtlapse, "zero")) or
      not (CCTK_EQUALS (initial_dtshift, "Kadath") or
           CCTK_EQUALS (initial_dtshift, "zero") or
           CCTK_EQUALS (initial_dtshift, "none")))
  {
    CCTK_PARAMWARN ("The parameters ADMBase::initial_data, ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse, and ADMBase::initial_dtshift must all be set to the value \"Kadath\" or \"none\"");
  }
}
}