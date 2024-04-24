#include <algorithm>
#include <cctk.h>
#include <cctk_Arguments.h>

#include "../test.hxx"
#include "reconstruct.hxx"
#include "reconx_utils.hxx"

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

void AsterXTests::test_ppm(std::mt19937_64 &engine, int repetitions) {
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

  std::array<CCTK_REAL, num_cells> press{};
  std::array<CCTK_REAL, num_cells> vel_dir{};

  for (std::size_t i = 0; i < num_cells; i++) {
    cell_centers[i] = x0 + dx * (i + 0.5);
  }

  for (std::size_t i = 0; i < num_faces; i++) {
    cell_faces[i] = x0 + i * dx;
  }

  // TODO: Set appropriate velocities and pressures
  for (auto &vel : vel_dir) {
    vel = 0.0;
  }

  for (auto &p : press) {
    p = 0.0;
  }

  uniform_real_distribution<CCTK_REAL> real_distrib{-1.0, 1.0};

  // Values obtained from parameter file defaults
  const reconstruct_params_t reconstruct_params{
      .ppm_shock_detection = true,
      .ppm_zone_flattening = false,
      .poly_k = 100.0,
      .poly_gamma = 2.0,
      .ppm_eta1 = 20.0,
      .ppm_eta2 = 0.05,
      .ppm_eps = 0.33,
      .ppm_eps_shock = 0.01,
      .ppm_small = 1.0e-7,
      .ppm_omega1 = 0.75,
      .ppm_omega2 = 10.0,
  };

  for (int i = 0; i < repetitions; i++) {
    CCTK_VINFO("Testing PPM reconstruction of a constant line, rep. %i", i);
    {
      const auto yval{real_distrib(engine)};
      const auto gf_faces{constant_line(yval, cell_faces)};
      const auto gf_centers{constant_line(yval, cell_centers)};

      const auto [recon_left, recon_right] = ppm_reconstruct(
          gf_centers[0], gf_centers[1], gf_centers[2], gf_centers[4],
          gf_centers[5], gf_centers[6], press[0], press[1], press[2], press[4],
          press[5], press[6], vel_dir[1], vel_dir[2], vel_dir[4], vel_dir[5],
          false, reconstruct_params);

      const auto real_left{gf_faces[3]}, real_right{gf_faces[4]};

      CCTK_VINFO("Real: %f %f | Recon: %f %f", real_left, real_right,
                 recon_left, recon_right);
    }
  }
}