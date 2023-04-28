#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "../test.hxx"
#include "../utils.hxx"

#include <random>

void AsterXTests::test_corss_product() {
  using namespace Arith;
  using namespace AsterX;

  using std::mt19937_64;
  using std::uniform_real_distribution;

  // Construct symmetric matrix of random real elements
  mt19937_64 engine{random_seed};
  uniform_real_distribution<CCTK_REAL> distrib{-1.0, 1.0};

  const vec<CCTK_REAL, 3> a{distrib(engine), distrib(engine), distrib(engine)};

  const vec<CCTK_REAL, 3> b{distrib(engine), distrib(engine), distrib(engine)};

  CCTK_VINFO("Testing cross product");
  {
    bool success{true};
    const auto actual_cross{calc_cross_product(a, b)};
    const auto expected_cross_0{a(1) * b(2) - a(2) * b(1)};
    const auto expected_cross_1{a(2) * b(0) - a(0) * b(2)};
    const auto expected_cross_2{a(0) * b(1) - a(1) * b(0)};

    if (!isapprox(actual_cross(0), expected_cross_0)) {
      CCTK_VINFO("  \033[31;1mFAILED\033[0m. Reason: Cross product element 0 "
                 "expected to be %.16f, but instead got %.16f",
                 expected_cross_0, actual_cross(0));
      success &= false;
    }

    if (!isapprox(actual_cross(1), expected_cross_1)) {
      CCTK_VINFO("  \033[31;1mFAILED\033[0m. Reason: Cross product element 1 "
                 "expected to be %.16f, but instead got %.16f",
                 expected_cross_1, actual_cross(1));
      success &= false;
    }

    if (!isapprox(actual_cross(2), expected_cross_2)) {
      CCTK_VINFO("  \033[31;1mFAILED\033[0m. Reason: Cross product element 2 "
                 "expected to be %.16f, but instead got %.16f",
                 expected_cross_2, actual_cross(2));
      success &= false;
    }

    if (success) {
      CCTK_VINFO("  \033[32;1mPASSED\033[0m.");
    }
  }
}