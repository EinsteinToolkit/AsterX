#ifndef ASTERX_TESTS_HXX
#define ASTERX_TESTS_HXX

#include <cmath>
#include <limits>
#include <random>

namespace AsterXTests {

/**
 * Seed value of the random number generation engines.
 */
constexpr const int random_seed = 100;

/**
 * Test two floating point values for approximate equality.
 *
 * This function is based on
 * https://docs.julialang.org/en/v1/base/math/#Base.isapprox
 *
 * @param x The first number to compare.
 * @param y The second number to compare.
 * @param atol The absolute tolerance for the comparison.
 * @param rtol The relative tolerance of the comparison.
 * @return True if x ~ y, false otherwise.
 */
template <typename fp_type>
bool isapprox(fp_type x, fp_type y, fp_type atol = 0.0) {
  using std::abs;
  using std::max;
  using std::sqrt;

  const fp_type rtol =
      atol > 0.0 ? 0.0 : sqrt(std::numeric_limits<fp_type>::epsilon());
  return abs(x - y) <= max(atol, rtol * max(abs(x), abs(y)));
}

void test_contraction_smat_upvec(std::mt19937_64 &engine, int repetitions);
void test_contraction_upvec_downvec(std::mt19937_64 &engine, int repetitions);

void test_cross_product(std::mt19937_64 &engine, int repetitions);

void test_wlorentz(std::mt19937_64 &engine, int repetitions);

void test_minmod(std::mt19937_64 &engine, int repetitions);

void test_mp5(std::mt19937_64 &engine, int repetitions);

} // namespace AsterXTests

#endif // ASTERX_TESTS_HXX