#include <cctk.h>
#include <cctk_Arguments.h>

#include <AMReX.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuMemory.H>

#include <array>
#include <cmath>

#include "setup_eos.hxx"

namespace EOSX {
namespace {

void check(const char *name, const CCTK_REAL actual,
           const CCTK_REAL expected) {
  // Compare small thermal quantities relatively, and never accept NaNs.
  if (!std::isfinite(actual) || !std::isfinite(expected) ||
      fabs(actual - expected) > 2.0e-11 * fmax(1.0e-12, fabs(expected)))
    CCTK_VERROR("EOSX test %s: actual=%.16e expected=%.16e",
               name, actual, expected);
}

void test_table() {
  // A log-linear table with known answers and rho/Ye-dependent energy
  // bounds. Exercise the real interpolator and inverse without a data file.
  std::array<CCTK_REAL, 3> lr{log(1.0e-6), log(1.0e-4), log(1.0e-2)};
  std::array<CCTK_REAL, 3> lt{log(1.0e-3), log(1.0e-2), log(1.0e-1)};
  std::array<CCTK_REAL, 2> ye{0.1, 0.5};
  std::array<CCTK_REAL, 18 * NTABLES> data{};
  for (int k = 0; k < 2; ++k)
    for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 3; ++i) {
        const int offset = NTABLES * (i + 3 * (j + 3 * k));
        data[offset + eos_3p_tabulated3d::PRESS] = lr[i] + lt[j];
        data[offset + eos_3p_tabulated3d::EPS] =
            lt[j] + 0.1 * lr[i] + 0.2 * ye[k];
        data[offset + eos_3p_tabulated3d::S] =
            lt[j] - 0.1 * lr[i] + 0.2 * ye[k];
      }
  linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> interp(
      data.data(), {3, 3, 2}, lr.data(), lt.data(), ye.data());

  for (CCTK_REAL shift : {0.0, 0.05}) {
    eos_3p_tabulated3d eos;
    eos.interptable = &interp;
    eos.energy_shift = &shift;
    eos.rgrho = {exp(lr.front()), exp(lr.back())};
    eos.rgtemp = {exp(lt.front()), exp(lt.back())};
    eos.rgye = {ye.front(), ye.back()};
    eos.rgeps = eos.compute_eps_range_full_table();

    const auto energy = [=](CCTK_REAL rho, CCTK_REAL temp, CCTK_REAL Ye) {
      return temp * pow(rho, 0.1) * exp(0.2 * Ye) - shift;
    };
    check("global eps minimum", eos.rgeps.min,
          energy(eos.rgrho.min, eos.rgtemp.min, eos.rgye.min));
    check("global eps maximum", eos.rgeps.max,
          energy(eos.rgrho.max, eos.rgtemp.max, eos.rgye.max));
    if (shift > 0.0 && !(eos.rgeps.min < 0.0))
      CCTK_ERROR("EOSX test: negative physical energy was discarded");

    for (const CCTK_REAL rho : {eos.rgrho.min, 3.0e-4, eos.rgrho.max})
      for (const CCTK_REAL Ye : {eos.rgye.min, 0.3, eos.rgye.max}) {
        const auto er = eos.range_eps_from_rho_ye(rho, Ye);
        check("local eps minimum", er.min, energy(rho, eos.rgtemp.min, Ye));
        check("local eps maximum", er.max, energy(rho, eos.rgtemp.max, Ye));
        for (const CCTK_REAL temp : {eos.rgtemp.min, 0.007, eos.rgtemp.max}) {
          CCTK_REAL eps = eos.eps_from_rho_temp_ye(rho, temp, Ye);
          check("physical eps", eps, energy(rho, temp, Ye));
          check("T -> eps -> T", eos.temp_from_rho_eps_ye(rho, eps, Ye), temp);
          check("inverse physical eps", eps, energy(rho, temp, Ye));
          check("inverse pressure", eos.press_from_rho_eps_ye(rho, eps, Ye),
                rho * temp);
          const CCTK_REAL kappa = log(temp) - 0.1 * log(rho) + 0.2 * Ye;
          check("table kappa from T",
                eos.kappa_from_rho_temp_ye(rho, temp, Ye), kappa);
          check("table kappa from eps",
                eos.kappa_from_rho_eps_ye(rho, eps, Ye), kappa);
        }
        for (const CCTK_REAL input : {eos.rgeps.min - 1.0, er.min,
                                      er.max, eos.rgeps.max + 1.0}) {
          CCTK_REAL eps = input;
          const CCTK_REAL temp = eos.temp_from_rho_eps_ye(rho, eps, Ye);
          const CCTK_REAL expected = fmin(fmax(input, er.min), er.max);
          check("inverse energy bound", eps, expected);
          check("eps -> T -> eps", eos.eps_from_rho_temp_ye(rho, temp, Ye),
                expected);
          check("inverse temperature bound", temp,
                input <= er.min ? eos.rgtemp.min : eos.rgtemp.max);
        }
      }
  }
}

void test_beq() {
  std::array<CCTK_REAL, 3> lr{log(1.0e-6), log(1.0e-4), log(1.0e-2)};
  std::array<CCTK_REAL, 3> lt{log(1.0e-3), log(1.0e-2), log(1.0e-1)};
  std::array<CCTK_REAL, 5> ye{0.1, 0.2, 0.3, 0.4, 0.5};
  std::array<CCTK_REAL, 45 * NTABLES> data{};
  linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> interp(
      data.data(), {3, 3, 5}, lr.data(), lt.data(), ye.data());
  eos_3p_tabulated3d eos;
  eos.interptable = &interp;
  eos.rgrho = {exp(lr.front()), exp(lr.back())};
  eos.rgtemp = {exp(lt.front()), exp(lt.back())};
  eos.rgye = {ye.front(), ye.back()};

  const auto fill = [&](const auto &func) {
    for (int k = 0; k < 5; ++k)
      for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i) {
          const int n = NTABLES * (i + 3 * (j + 3 * k));
          data[n + eos_3p_tabulated3d::MU_E] =
              1.0 + func(lr[i], lt[j], ye[k]);
          data[n + eos_3p_tabulated3d::MU_P] = 1.0;
          data[n + eos_3p_tabulated3d::MU_N] = 2.0;
        }
  };

  // Exercise off-grid roots, both bracket orientations and bounded rho/T.
  for (CCTK_REAL sign : {-1.0, 1.0}) {
    fill([&](CCTK_REAL r, CCTK_REAL t, CCTK_REAL y) {
      return sign * (y - (0.3 + 0.01 * (r - lr[1]) +
                          0.02 * (t - lt[1])));
    });
    for (CCTK_REAL rho : {0.5 * eos.rgrho.min, 3.0e-4,
                          2.0 * eos.rgrho.max})
      for (CCTK_REAL temp : {0.5 * eos.rgtemp.min, 0.007,
                             2.0 * eos.rgtemp.max}) {
        const CCTK_REAL r = std::clamp(rho, eos.rgrho.min, eos.rgrho.max);
        const CCTK_REAL t = std::clamp(temp, eos.rgtemp.min, eos.rgtemp.max);
        check("beta-equilibrium Ye", eos.ye_beq_from_rho_temp(rho, temp),
              0.3 + 0.01 * (log(r) - lr[1]) +
                  0.02 * (log(t) - lt[1]));
      }
  }

  // Exact roots and the documented nearest-endpoint fallback.
  for (CCTK_REAL target : {0.1, 0.2, 0.5, 0.0, 0.8}) {
    fill([&](CCTK_REAL, CCTK_REAL, CCTK_REAL y) { return y - target; });
    const CCTK_REAL expected =
        std::clamp(target, eos.rgye.min, eos.rgye.max);
    check("bounded beta-equilibrium Ye",
          eos.ye_beq_from_rho_temp(1.0e-4, 0.01), expected);
  }

  // Use the actual piecewise-linear table, not one line across all Ye points.
  fill([](CCTK_REAL, CCTK_REAL, CCTK_REAL y) { return y * y - 0.13; });
  check("piecewise beta-equilibrium Ye",
        eos.ye_beq_from_rho_temp(1.0e-4, 0.01),
        0.3 + 0.1 * (0.13 - 0.09) / (0.16 - 0.09));
}

void test_ideal() {
  for (const CCTK_REAL gamma : {1.4, 2.0}) {
    eos_3p_idealgas eos;
    eos_3p::range er{0.0, 2.0}, rr{1.0e-12, 1.0}, yr{0.0, 1.0};
    eos.init(gamma, 2.0, er, rr, yr);
    for (const CCTK_REAL rho : {1.0e-10, 1.0e-4, 0.5})
      for (CCTK_REAL eps : {0.0, 1.0e-8, 0.1, 2.0}) {
        const CCTK_REAL temp = eos.temp_from_rho_eps_ye(rho, eps, 0.5);
        check("ideal temperature", temp, 2.0 * (gamma - 1.0) * eps);
        check("ideal energy", eos.eps_from_rho_temp_ye(rho, temp, 0.5), eps);
        check("ideal pressure", eos.press_from_rho_eps_ye(rho, eps, 0.5),
              (gamma - 1.0) * rho * eps);
        const CCTK_REAL kappa =
            (gamma - 1.0) * eps * pow(rho, 1.0 - gamma);
        check("ideal kappa from T",
              eos.kappa_from_rho_temp_ye(rho, temp, 0.5), kappa);
        check("ideal kappa from eps",
              eos.kappa_from_rho_eps_ye(rho, eps, 0.5), kappa);
      }
  }
}

template <typename EOSType> void test_device(const EOSType *eos) {
  // Only the active EOS has device-accessible storage. Synthetic tables
  // above remain host-local; never capture their pointers in a device loop.
  amrex::Gpu::DeviceScalar<unsigned int> failures(0);
  auto *failed = failures.dataPtr();
  amrex::ParallelFor(32, [=] AMREX_GPU_DEVICE(int i) {
    const CCTK_REAL f = (i + 0.5) / 32.0;
    const CCTK_REAL rho =
        exp((1.0 - f) * log(eos->rgrho.min) + f * log(eos->rgrho.max));
    const CCTK_REAL temp =
        eos->rgtemp.min + f * (eos->rgtemp.max - eos->rgtemp.min);
    const CCTK_REAL Ye = eos->rgye.min + f * (eos->rgye.max - eos->rgye.min);
    const CCTK_REAL eps = eos->eps_from_rho_temp_ye(rho, temp, Ye);
    CCTK_REAL eps_back = eps;
    const CCTK_REAL temp_back = eos->temp_from_rho_eps_ye(rho, eps_back, Ye);
    const CCTK_REAL press = eos->press_from_rho_temp_ye(rho, temp, Ye);
    const CCTK_REAL press_back = eos->press_from_rho_eps_ye(rho, eps_back, Ye);
    const CCTK_REAL kappa = eos->kappa_from_rho_temp_ye(rho, temp, Ye);
    const CCTK_REAL kappa_back = eos->kappa_from_rho_eps_ye(rho, eps_back, Ye);
    if (!std::isfinite(eps) || !std::isfinite(eps_back) ||
        !std::isfinite(temp_back) || !std::isfinite(press) ||
        !std::isfinite(press_back) || !std::isfinite(kappa) ||
        !std::isfinite(kappa_back) ||
        fabs(temp_back - temp) > 1.0e-7 * fmax(temp, 1.0e-12) ||
        fabs(eps_back - eps) > 1.0e-7 * fmax(fabs(eps), 1.0e-12) ||
        fabs(press_back - press) > 1.0e-7 * fmax(fabs(press), 1.0e-20) ||
        fabs(kappa_back - kappa) > 1.0e-7 * fmax(fabs(kappa), 1.0))
      amrex::HostDevice::Atomic::Add(failed, 1U);
  });
  amrex::Gpu::streamSynchronize();
  if (failures.dataValue())
    CCTK_ERROR("EOSX test: active-EOS device round trip failed");
}
} // namespace

extern "C" void EOSX_Test(CCTK_ARGUMENTS) {
  test_table();
  test_beq();
  test_ideal();
  if (global_eos_3p_ig)
    test_device(global_eos_3p_ig);
  if (global_eos_3p_tab3d)
    test_device(global_eos_3p_tab3d);
  CCTK_INFO("EOSX energy-bound, temperature and beta-equilibrium tests passed");
}

} // namespace EOSX
