#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "../test.hxx"
#include "aster_utils.hxx"

#include <random>

void AsterXTests::test_contraction_upvec_downvec(std::mt19937_64 &engine,
                                                 int repetitions) {
  using namespace Arith;
  using namespace AsterUtils;
  using std::uniform_real_distribution;

  uniform_real_distribution<CCTK_REAL> distrib{-1.0, 1.0};

  for (int i = 0; i < repetitions; i++) {

    const vec<CCTK_REAL, 3> v_up{distrib(engine), distrib(engine),
                                 distrib(engine)};
    const vec<CCTK_REAL, 3> v_down{distrib(engine), distrib(engine),
                                   distrib(engine)};

    CCTK_VINFO("Testing up_vec - down_vec contraction repetition %i", i);
    {
      bool success{true};
      const auto actual_contraction{calc_contraction(v_up, v_down)};
      const auto expected_contraction{
          v_up(0) * v_down(0) + v_up(1) * v_down(1) + v_up(2) * v_down(2)};

      if (!isapprox(actual_contraction, expected_contraction)) {
        CCTK_VINFO("  FAILED. Reason: Contraction expected to "
                   "be %.16f, but instead got %.16f",
                   expected_contraction, actual_contraction);
        success &= false;
      }

      if (success) {
        CCTK_VINFO("  PASSED.");
      }
    }
  }
}
