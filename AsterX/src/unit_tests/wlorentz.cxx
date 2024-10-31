#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "../test.hxx"
#include "aster_utils.hxx"

// Value must be bounded
void AsterXTests::test_wlorentz(std::mt19937_64 &engine, int repetitions) {
  using namespace Arith;
  using namespace AsterUtils;
  using std::sqrt;
  using std::uniform_real_distribution;

  uniform_real_distribution<CCTK_REAL> distrib{0.0, 0.9};

  for (int i = 0; i < repetitions; i++) {

    const vec<CCTK_REAL, 3> v_up{distrib(engine), distrib(engine),
                                 distrib(engine)};

    const vec<CCTK_REAL, 3> v_down{distrib(engine), distrib(engine),
                                   distrib(engine)};

    const auto contraction{calc_contraction(v_up, v_down)};

    CCTK_VINFO("Testing wlorentz repetition %i", i);
    {
      bool success{true};

      const auto actual_lorentz{calc_wlorentz(v_up, v_down)};
      const auto expected_lorentz{1.0 / sqrt(1.0 - contraction)};

      if (!isapprox(actual_lorentz, expected_lorentz)) {
        CCTK_VINFO("  FAILED. Reason: expected %.f16 and got %.f16",
                   expected_lorentz, actual_lorentz);
        success &= false;
      }

      if (success) {
        CCTK_VINFO("  PASSED.");
      }
    }
  }
}
