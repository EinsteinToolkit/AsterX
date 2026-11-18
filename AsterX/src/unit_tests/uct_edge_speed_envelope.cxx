#include <cctk.h>

#include "../fluxes.hxx"
#include "../test.hxx"

void AsterXTests::test_uct_edge_speed_envelope() {
  using namespace AsterX;

  // Each adjacent face supplies one of the two extrema. Either one-face
  // sampling would fail this case.
  const auto speeds = uct_edge_speed_envelope(
      CCTK_REAL(0.25), CCTK_REAL(4.0), CCTK_REAL(3.0), CCTK_REAL(0.5));
  if (speeds(0) != CCTK_REAL(3.0) || speeds(1) != CCTK_REAL(4.0))
    CCTK_ERROR("UCT edge speed envelope omitted an adjacent face");

  // Face ordering must not affect an edge-centered value.
  const auto reversed = uct_edge_speed_envelope(
      CCTK_REAL(3.0), CCTK_REAL(0.5), CCTK_REAL(0.25), CCTK_REAL(4.0));
  if (reversed(0) != speeds(0) || reversed(1) != speeds(1))
    CCTK_ERROR("UCT edge speed envelope depends on face ordering");

  CCTK_INFO("UCT two-face edge speed envelope passed");
}

void AsterXTests::test_hlle() {
  using namespace Arith;
  using namespace AsterX;

  const vec<CCTK_REAL, 4> zero{0.0, 0.0, 0.0, 0.0};
  const vec<vec<CCTK_REAL, 4>, 2> zero_speeds{zero, zero};
  const vec<CCTK_REAL, 2> var{1.0, 3.0};
  const vec<CCTK_REAL, 2> flux{2.0, 4.0};

  if (hlle(zero_speeds, var, flux) != CCTK_REAL(3.0))
    CCTK_ERROR("HLLE degenerate-speed fallback failed");

  const vec<vec<CCTK_REAL, 4>, 2> speeds{
      vec<CCTK_REAL, 4>{-2.0, -1.0, 1.0, 2.0},
      vec<CCTK_REAL, 4>{-1.0, 0.0, 2.0, 3.0}};
  if (!isapprox(hlle(speeds, var, flux), CCTK_REAL(0.4)))
    CCTK_ERROR("HLLE flux calculation failed");

  CCTK_INFO("HLLE tests passed");
}
