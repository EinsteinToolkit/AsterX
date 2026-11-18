#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include <array>
#include <cmath>

#include "seeds_utils.hxx"
#include "setup_eos.hxx"

namespace AsterSeeds {
using namespace std;
using namespace EOSX;

extern "C" void AsterSeeds_TestBetaFloor(CCTK_ARGUMENTS) {
  std::array<CCTK_REAL, 2> lr{log(1.0e-6), log(1.0e-2)};
  std::array<CCTK_REAL, 2> lt{log(1.0e-3), log(1.0e-1)};
  std::array<CCTK_REAL, 2> ye{0.1, 0.5};
  std::array<CCTK_REAL, 8 * NTABLES> data{};
  for (int k = 0; k < 2; ++k)
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < 2; ++i) {
        const int n = NTABLES * (i + 2 * (j + 2 * k));
        data[n + eos_3p_tabulated3d::PRESS] = lr[i] + lt[j];
        data[n + eos_3p_tabulated3d::EPS] = lt[j];
        data[n + eos_3p_tabulated3d::S] = lt[j] - lr[i];
      }
  linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> interp(
      data.data(), {2, 2, 2}, lr.data(), lt.data(), ye.data());
  CCTK_REAL shift = 0.05;
  eos_3p_tabulated3d eos;
  eos.interptable = &interp;
  eos.energy_shift = &shift;
  eos.rgrho = {exp(lr[0]), exp(lr[1])};
  eos.rgtemp = {exp(lt[0]), exp(lt[1])};
  eos.rgye = {ye[0], ye[1]};
  eos.rgeps = eos.compute_eps_range_full_table();
  CCTK_REAL rho = 0.0, eps = 0.0, press = 0.0, entropy = 0.0;
  if (!set_beta_floor(&eos, 6.0e-6, 0.02, 0.3, rho, eps, press, entropy) ||
      fabs(rho - 3.0e-4) > 1.0e-14 || fabs(eps + 0.03) > 1.0e-12 ||
      fabs(press - 6.0e-6) > 1.0e-15 ||
      fabs(entropy - log(0.02 / 3.0e-4)) > 1.0e-10)
    CCTK_ERROR("AsterSeeds test: inconsistent beta-floor state");
  const CCTK_REAL rho0 = rho, eps0 = eps, press0 = press, entropy0 = entropy;
  if (set_beta_floor(&eos, 1.0, 0.02, 0.3, rho, eps, press, entropy) ||
      rho != rho0 || eps != eps0 || press != press0 || entropy != entropy0)
    CCTK_ERROR("AsterSeeds test: unattainable beta floor changed the state");
  // A steep pressure interpolant amplifies the density-inversion roundoff.
  std::array<CCTK_REAL, 3> lr_small{log(1.0e-20), log(1.0e-12), log(1.0e-4)};
  std::array<CCTK_REAL, 12 * NTABLES> small_data{};
  for (int k = 0; k < 2; ++k)
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < 3; ++i) {
        const int n = NTABLES * (i + 3 * (j + 2 * k));
        small_data[n + eos_3p_tabulated3d::PRESS] = 3.0 * lr_small[i] + lt[j];
        small_data[n + eos_3p_tabulated3d::EPS] = lt[j];
      }
  linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> small(
      small_data.data(), {3, 2, 2}, lr_small.data(), lt.data(), ye.data());
  eos.interptable = &small;
  eos.rgrho = {exp(lr_small[0]), exp(lr_small[2])};
  const CCTK_REAL rhoL = 1.624501063509632e-14;
  const CCTK_REAL target = eos.press_from_rho_temp_ye(rhoL, 0.02, 0.3);
  if (!set_beta_floor(&eos, target, 0.02, 0.3, rho, eps, press, entropy) ||
      fabs(rho / rhoL - 1.0) > 1.0e-12 ||
      fabs(press / target - 1.0) > 1.0e-12)
    CCTK_ERROR("AsterSeeds test: attainable beta floor rejected");
  CCTK_INFO("AsterSeeds beta-floor tests passed");
}

} // namespace AsterSeeds
