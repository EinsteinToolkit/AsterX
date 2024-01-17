#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "test.hxx"

extern "C" void AsterX_Test(CCTK_ARGUMENTS) {
  using namespace AsterXTests;

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (unit_test) {
    std::mt19937_64 engine{random_seed};
    const int repetitions{unit_test_repetitions};

    CCTK_VINFO("Running unit tests with %i repetitions", repetitions);

    test_contraction_smat_upvec(engine, repetitions);
    test_contraction_upvec_downvec(engine, repetitions);

    test_cross_product(engine, repetitions);

    test_wlorentz(engine, repetitions);

    test_minmod(engine, repetitions);

    test_mp5(engine, repetitions);

  } else {
    CCTK_INFO("Skipping unit tests");
  }
}
