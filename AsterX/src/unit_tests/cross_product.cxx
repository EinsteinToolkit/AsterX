#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "../test.hxx"
#include "aster_utils.hxx"

#include <random>

void AsterXTests::test_cross_product(std::mt19937_64 &engine, int repetitions) {
  using namespace Arith;
  using namespace AsterUtils;
  using std::uniform_real_distribution;

  uniform_real_distribution<CCTK_REAL> distrib{-1.0, 1.0};

  for (int i = 0; i < repetitions; i++) {

    const vec<CCTK_REAL, 3> a{distrib(engine), distrib(engine),
                              distrib(engine)};
    const vec<CCTK_REAL, 3> b{distrib(engine), distrib(engine),
                              distrib(engine)};

    CCTK_VINFO("Testing cross product repetition %i", i);
    {
      bool success{true};
      const auto actual_cross{calc_cross_product(a, b)};
      const auto expected_cross_0{a(1) * b(2) - a(2) * b(1)};
      const auto expected_cross_1{a(2) * b(0) - a(0) * b(2)};
      const auto expected_cross_2{a(0) * b(1) - a(1) * b(0)};

      if (!isapprox(actual_cross(0), expected_cross_0)) {
        CCTK_VINFO("  FAILED. Reason: Cross product element 0 "
                   "expected to be %.16f, but instead got %.16f",
                   expected_cross_0, actual_cross(0));
        success &= false;
      }

      if (!isapprox(actual_cross(1), expected_cross_1)) {
        CCTK_VINFO("  FAILED. Reason: Cross product element 1 "
                   "expected to be %.16f, but instead got %.16f",
                   expected_cross_1, actual_cross(1));
        success &= false;
      }

      if (!isapprox(actual_cross(2), expected_cross_2)) {
        CCTK_VINFO("  FAILED. Reason: Cross product element 2 "
                   "expected to be %.16f, but instead got %.16f",
                   expected_cross_2, actual_cross(2));
        success &= false;
      }

      if (success) {
        CCTK_VINFO("  PASSED.");
      }
    }
  }
}
