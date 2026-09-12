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
