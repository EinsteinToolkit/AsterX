#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "../test.hxx"
#include "aster_utils.hxx"

#include <random>

void AsterXTests::test_contraction_smat_upvec(std::mt19937_64 &engine,
                                              int repetitions) {
  using namespace Arith;
  using namespace AsterUtils;
  using std::uniform_real_distribution;

  uniform_real_distribution<CCTK_REAL> distrib{-1.0, 1.0};

  for (int i = 0; i < repetitions; i++) {

    const smat<CCTK_REAL, 3> m{distrib(engine), distrib(engine),
                               distrib(engine), distrib(engine),
                               distrib(engine), distrib(engine)};
    const vec<CCTK_REAL, 3> v_up{distrib(engine), distrib(engine),
                                 distrib(engine)};

    CCTK_VINFO("Testing smat - up_vec contraction repetition %i", i);
    {
      bool success{true};
      const auto actual_contraction{calc_contraction(m, v_up)};
      const auto expected_contraction_0{m(0, 0) * v_up(0) + m(0, 1) * v_up(1) +
                                        m(0, 2) * v_up(2)};
      const auto expected_contraction_1{m(1, 0) * v_up(0) + m(1, 1) * v_up(1) +
                                        m(1, 2) * v_up(2)};
      const auto expected_contraction_2{m(2, 0) * v_up(0) + m(2, 1) * v_up(1) +
                                        m(2, 2) * v_up(2)};

      if (!isapprox(actual_contraction(0), expected_contraction_0)) {
        CCTK_VINFO("  FAILED. Reason: Contraction component 0 "
                   "expected to be %.16f, but instead got %.16f",
                   expected_contraction_0, actual_contraction(0));
        success &= false;
      }

      if (!isapprox(actual_contraction(1), expected_contraction_1)) {
        CCTK_VINFO("  FAILED. Reason: Contraction component 1 "
                   "expected to be %.16f, but instead got %.16f",
                   expected_contraction_0, actual_contraction(0));
        success &= false;
      }

      if (!isapprox(actual_contraction(2), expected_contraction_2)) {
        CCTK_VINFO("  FAILED. Reason: Contraction component 2 "
                   "expected to be %.16f, but instead got %.16f",
                   expected_contraction_0, actual_contraction(0));
        success &= false;
      }

      if (success) {
        CCTK_VINFO("  PASSED.");
      }
    }
  }
}
