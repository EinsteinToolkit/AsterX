#include <cctk.h>
#include <cctk_Arguments.h>

#include "../test.hxx"
#include "reconstruct.hxx"

#include <random>

void AsterXTests::test_minmod(std::mt19937_64 &engine, int repetitions) {
  using namespace ReconX;
  using std::abs;
  using std::max;
  using std::min;
  using std::uniform_real_distribution;

  uniform_real_distribution<CCTK_REAL> positive_distrib{0.0, 1.0};
  uniform_real_distribution<CCTK_REAL> negative_distrib{-1.0, 0.0};
  uniform_real_distribution<CCTK_REAL> real_distrib{-1.0, 1.0};

  for (int i = 0; i < repetitions; i++) {

    const auto p1{positive_distrib(engine)};
    const auto p2{positive_distrib(engine)};
    const auto p3{positive_distrib(engine)};
    const auto p4{positive_distrib(engine)};

    const auto n1{negative_distrib(engine)};
    const auto n2{negative_distrib(engine)};
    const auto n3{negative_distrib(engine)};
    const auto n4{negative_distrib(engine)};

    const auto r1{real_distrib(engine)};
    const auto r2{real_distrib(engine)};

    const auto positive_minmod{minmod(p1, p2)};
    const auto negative_minmod{minmod(n1, n2)};

    const auto mixed_minmod_1{minmod(n1, p1)};
    const auto mixed_minmod_2{minmod(n2, p2)};

    const auto real_minmod_1{minmod(r1, r2)};

    const auto min_p1p2{min(p1, p2)};
    const auto max_n1n2{max(n1, n2)};

    const auto max_abs{max(abs(r1), abs(r2))};

    CCTK_VINFO("Testing minmod of positive numbers, rep. %i", i);
    {

      if (!isapprox(positive_minmod, min_p1p2)) {
        CCTK_VINFO("  FAILED. Reason: minmod with positive "
                   "numbers %.f16 and %.f16 expected to be %.f16 but resulted "
                   "in %.f16",
                   p1, p2, min_p1p2, positive_minmod);
      } else {
        CCTK_VINFO("  PASSED.");
      }
    }

    CCTK_VINFO("Testing minmod of negative numbers, rep. %i", i);
    {
      if (!isapprox(negative_minmod, max_n1n2)) {
        CCTK_VINFO("  FAILED. Reason: minmod with negative "
                   "numbers %.f16 and %.f16 expected to be %.f16 but resulted "
                   "in %.f16",
                   n1, n2, max_n1n2, negative_minmod);
      } else {
        CCTK_VINFO("  PASSED.");
      }
    }

    CCTK_VINFO("Testing minmod of mixed sign numbers, rep. %i", i);
    {
      if (!isapprox(mixed_minmod_1, CCTK_REAL(0))) {
        CCTK_VINFO(
            "  FAILED. Reason: minmod with mixed sign "
            "numbers %.f16 and %.f16 expected tp be 0.0 but resulted in %.f16",
            n1, p1, mixed_minmod_1);
      } else {
        CCTK_VINFO("  PASSED.");
      }

      if (!isapprox(mixed_minmod_1, CCTK_REAL(0))) {
        CCTK_VINFO(
            "  FAILED. Reason: minmod with mixed sign "
            "numbers %.f16 and %.f16 expected tp be 0.0 but resulted in %.f16",
            n2, p2, mixed_minmod_2);
      } else {
        CCTK_VINFO("  PASSED.");
      }
    }

    CCTK_VINFO("Testing minmod(a, b) <= max(|a|,|b|), rep. %i", i);
    {
      if (!(real_minmod_1 < max_abs || isapprox(real_minmod_1, max_abs))) {
        CCTK_VINFO("  FAILED. Reason: minmod with reals %.f16 and %.f16 "
                   "returned %.f, which is not bounded by %.f16",
                   r1, r2, real_minmod_1, max_abs);
      } else {
        CCTK_VINFO("  PASSED.");
      }
    }

    CCTK_VINFO("Testing minmod(a, b, c, d) = minmod(minmod(a, b), minmod(c, "
               "d)) for positive arguments, rep. %i",
               i);
    {
      const auto expected_min{min(p1, min(p2, min(p3, p4)))};
      const auto computed_min(minmod(minmod(p1, p2), minmod(p3, p4)));

      if (!isapprox(expected_min, computed_min)) {
        CCTK_VINFO("  FAILED. Reason: expected %.f16 and got %.f16",
                   expected_min, computed_min);
      } else {
        CCTK_VINFO("  PASSED.");
      }
    }

    CCTK_VINFO("Testing minmod(a, b, c, d) = minmod(minmod(a, b), minmod(c, "
               "d)) for negative arguments, rep. %i",
               i);
    {
      const auto expected_max{max(n1, max(n2, max(n3, n4)))};
      const auto computed_max(minmod(minmod(n1, n2), minmod(n3, n4)));

      if (!isapprox(expected_max, computed_max)) {
        CCTK_VINFO("  FAILED. Reason: expected %.f16 and got %.f16",
                   expected_max, computed_max);
      } else {
        CCTK_VINFO("  PASSED.");
      }
    }
  }
}