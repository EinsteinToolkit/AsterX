#include <loop_device.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

// initial setup adapted from Springel+2010

namespace KHInitial {

extern "C" void KHInitial_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_KHInitial_Initialize;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_hydro, "KHI")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          // set velocities in opposite directions in x along the slip surface
          // set slightly different densities along the slip surface
          // location of slip surfaces hardcoded at (abs(p.y) = 0.25)
          using std::abs;
          if (abs(p.y) >= 0.25) {
            rho(p.I) = rhoUp;
            velx(p.I) = vxUp;
          } else {
            rho(p.I) = rhoLow;
            velx(p.I) = vxLow;
          }
          // excite the instability by peturbing v^y
          using std::exp, std::pow, std::sin;
          vely(p.I) = w0 * sin(4 * M_PI * p.x) *
                      (exp(-pow(p.y - 0.25, 2) / (2 * pow(sigma, 2))) +
                       exp(-pow(p.y + 0.25, 2) / (2 * pow(sigma, 2))));
          velz(p.I) = 0.0;

          // set constant initial pressure throughout the domain
          press(p.I) = p_val;

          // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          eps(p.I) = press(p.I) / (rho(p.I) * (gamma - 1));
        });

  } else {
    CCTK_ERROR("Internal error");
  }
}

} // namespace KHInitial
