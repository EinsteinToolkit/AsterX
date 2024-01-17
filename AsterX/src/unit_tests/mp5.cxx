#include <algorithm>
#include <cctk.h>
#include <cctk_Arguments.h>

#include "../test.hxx"
#include "reconstruct.hxx"

#include <array>

template <std::size_t N>
static std::array<CCTK_REAL, N>
constant_line(CCTK_REAL yvalue, const std::array<CCTK_REAL, N> &) {
  std::array<CCTK_REAL, N> func{};
  std::fill(func.begin(), func.end(), yvalue);
  return func;
}

template <std::size_t N>
static std::array<CCTK_REAL, N> sloped_line(CCTK_REAL a, CCTK_REAL b,
                                            const std::array<CCTK_REAL, N> &x) {
  std::array<CCTK_REAL, N> func{};

  for (std::size_t i = 0; i < N; i++) {
    func[i] = a * x[i] + b;
  }

  return func;
}

template <std::size_t N>
static std::array<CCTK_REAL, N> parabola(CCTK_REAL a, CCTK_REAL b, CCTK_REAL c,
                                         const std::array<CCTK_REAL, N> &x) {
  std::array<CCTK_REAL, N> func{};

  for (std::size_t i = 0; i < N; i++) {
    func[i] = a * x[i] * x[i] + b * x[i] + c;
  }

  return func;
}

template <std::size_t N>
static std::array<CCTK_REAL, N> hyperbole(CCTK_REAL a, CCTK_REAL b, CCTK_REAL c,
                                          CCTK_REAL d,
                                          const std::array<CCTK_REAL, N> &x) {
  std::array<CCTK_REAL, N> func{};

  for (std::size_t i = 0; i < N; i++) {
    func[i] = a * x[i] * x[i] * x[i] + b * x[i] * x[i] + c * x[i] + d;
  }

  return func;
}

void AsterXTests::test_mp5(std::mt19937_64 &engine, int repetitions) {
  using namespace ReconX;
  using std::uniform_real_distribution;

  // Test domain
  static constexpr const std::size_t num_cells{7};
  static constexpr const std::size_t num_faces{num_cells + 1};
  static constexpr const CCTK_REAL x0{-1.0};
  static constexpr const CCTK_REAL xf{1.0};
  static constexpr const CCTK_REAL dx{(xf - x0) / num_cells};

  // Cell faces and centers
  std::array<CCTK_REAL, num_faces> cell_faces{};
  std::array<CCTK_REAL, num_cells> cell_centers{};

  for (std::size_t i = 0; i < num_cells; i++) {
    cell_centers[i] = x0 + dx * (i + 0.5);
  }

  for (std::size_t i = 0; i < num_faces; i++) {
    cell_faces[i] = x0 + i * dx;
  }

  uniform_real_distribution<CCTK_REAL> real_distrib{-1.0, 1.0};

  for (int i = 0; i < repetitions; i++) {

    CCTK_VINFO("Testing mp5 reconstruction of a constant line, rep. %i", i);
    {
      const auto yval{real_distrib(engine)};
      const auto gf_faces{constant_line(yval, cell_faces)};
      const auto gf_centers{constant_line(yval, cell_centers)};

      const auto [recon_left, recon_right] =
          mp5_reconstruct(gf_centers[0], gf_centers[1], gf_centers[2],
                          gf_centers[4], gf_centers[5], gf_centers[6], 4.0);

      const auto real_left{gf_faces[3]}, real_right{gf_faces[4]};

      if (!isapprox(recon_left, real_left)) {
        CCTK_VINFO(
            "  FAILED. Reason: Reconstructed values %.16f are not "
            "comparable to the true values %.16f at the left "
            "face. This means that the reconstruction is introducing new "
            "global maxima or minima to the reconstructed function",
            recon_left, real_left);

      } else if (!isapprox(recon_right, real_right)) {
        CCTK_VINFO(
            "  FAILED. Reason: Reconstructed values %.16f are not "
            "comparable to the  values %.16f at the right "
            "face. This means that the reconstruction is introducing new "
            "global maxima or minima to the reconstructed function",
            recon_right, real_right);
      } else {
        CCTK_VINFO("  PASSED.");
      }
    }

    CCTK_VINFO("Testing mp5 reconstruction of a sloped line rep. %i", i);
    {
      const auto a{real_distrib(engine)};
      const auto b{real_distrib(engine)};
      const auto gf_faces{sloped_line(a, b, cell_faces)};
      const auto gf_centers{sloped_line(a, b, cell_centers)};

      const auto [recon_left, recon_right] =
          mp5_reconstruct(gf_centers[0], gf_centers[1], gf_centers[2],
                          gf_centers[4], gf_centers[5], gf_centers[6], 4.0);

      const auto global_max{
          *std::max_element(gf_faces.begin(), gf_faces.end())};

      const auto global_min{
          *std::min_element(gf_faces.begin(), gf_faces.end())};

      if (recon_left > global_max || recon_left < global_min) {
        CCTK_VINFO(
            "  FAILED. Reason: Left reconstructed value %.16f introduces new "
            "global minma (%.16f) or maxima (%.16f).",
            recon_left, global_min, global_max);
      } else if (recon_right > global_max || recon_right < global_min) {
        CCTK_VINFO(
            "  FAILED. Reason: Right reconstructed value %.16f introduces new "
            "global minma (%.16f) or maxima (%.16f).",
            recon_right, global_min, global_max);
      } else {
        CCTK_VINFO("  PASSED.");
      }
    }

    CCTK_VINFO("Testing mp5 reconstruction of a parabola rep. %i", i);
    {
      const auto a{real_distrib(engine)};
      const auto b{real_distrib(engine)};
      const auto c{real_distrib(engine)};
      const auto gf_faces{parabola(a, b, c, cell_faces)};
      const auto gf_centers{parabola(a, b, c, cell_centers)};

      const auto [recon_left, recon_right] =
          mp5_reconstruct(gf_centers[0], gf_centers[1], gf_centers[2],
                          gf_centers[4], gf_centers[5], gf_centers[6], 4.0);

      const auto global_max{
          *std::max_element(gf_faces.begin(), gf_faces.end())};

      const auto global_min{
          *std::min_element(gf_faces.begin(), gf_faces.end())};

      if (recon_left > global_max || recon_left < global_min) {
        CCTK_VINFO(
            "  FAILED. Reason: Left reconstructed value %.16f introduces new "
            "global minma (%.16f) or maxima (%.16f).",
            recon_left, global_min, global_max);
      } else if (recon_right > global_max || recon_right < global_min) {
        CCTK_VINFO(
            "  FAILED. Reason: Right reconstructed value %.16f introduces new "
            "global minma (%.16f) or maxima (%.16f).",
            recon_right, global_min, global_max);
      } else {
        CCTK_VINFO("  PASSED.");
      }
    }

    CCTK_VINFO("Testing mp5 reconstruction of a hyperbole rep. %i", i);
    {
      const auto a{real_distrib(engine)};
      const auto b{real_distrib(engine)};
      const auto c{real_distrib(engine)};
      const auto d{real_distrib(engine)};
      const auto gf_faces{hyperbole(a, b, c, d, cell_faces)};
      const auto gf_centers{hyperbole(a, b, c, d, cell_centers)};

      const auto [recon_left, recon_right] =
          mp5_reconstruct(gf_centers[0], gf_centers[1], gf_centers[2],
                          gf_centers[4], gf_centers[5], gf_centers[6], 4.0);

      const auto global_max{
          *std::max_element(gf_faces.begin(), gf_faces.end())};

      const auto global_min{
          *std::min_element(gf_faces.begin(), gf_faces.end())};

      if (recon_left > global_max || recon_left < global_min) {
        CCTK_VINFO(
            "  FAILED. Reason: Left reconstructed value %.16f introduces new "
            "global minma (%.16f) or maxima (%.16f).",
            recon_left, global_min, global_max);
      } else if (recon_right > global_max || recon_right < global_min) {
        CCTK_VINFO(
            "  FAILED. Reason: Right reconstructed value %.16f introduces new "
            "global minma (%.16f) or maxima (%.16f).",
            recon_right, global_min, global_max);
      } else {
        CCTK_VINFO("  PASSED.");
      }
    }
  }
}