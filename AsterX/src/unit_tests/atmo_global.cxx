#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuMemory.H>

#include <limits>

#include "../atmo_global.hxx"
#include "setup_eos.hxx"

namespace AsterX {
using namespace EOSX;
using namespace Con2PrimFactory;

namespace {

void CheckAtmoCopy(const atmosphere &a, const atmosphere &b) {
  if (a.rho_atmo != b.rho_atmo || a.eps_atmo != b.eps_atmo ||
      a.ye_atmo != b.ye_atmo || a.press_atmo != b.press_atmo ||
      a.temp_atmo != b.temp_atmo || a.entropy_atmo != b.entropy_atmo ||
      a.rho_cut != b.rho_cut)
    CCTK_ERROR("Global atmosphere copy changed the stored state");
}

void TestAtmoGlobal() {
  eos_1p_polytropic cold;
  cold.init(2.0, 100.0, 1.0);
  eos_1p_polytropic other_cold;
  other_cold.init(2.0, 100.0, 1.0);
  eos_3p_idealgas eos;
  eos_3p::range er{0.0, 1.0}, rr{1.0e-12, 1.0}, yr{0.0, 1.0};
  eos.init(2.0, 1.0, er, rr, yr);
  auto other_eos = eos;

  for (int mode = 0; mode < 3; ++mode) {
    bool thermal = mode != 0, use_press = mode == 2, Ye_beq = false;
    CCTK_REAL rho = 1.0e-6, press = 1.0e-10, temp = 0.02, Ye = 0.3;
    CCTK_REAL nr = 0.0, np = 0.0, nt = 0.0, tol = 0.001;
    const void *cold_ptr = &cold, *eos_ptr = &eos;
    atmo_global saved;
    atmosphere state{};
    const auto load = [&]() {
      return saved.load(cold_ptr, eos_ptr, rho, press, temp, Ye,
                        nr, np, nt, tol, thermal, use_press, Ye_beq, state);
    };
    if (load())
      CCTK_ERROR("Uninitialized global atmosphere was used");

    const auto direct =
        make_atmo(&cold, &eos, 0.0, rho, press, temp, Ye, 10.0,
                  nr, np, nt, tol, thermal, use_press, Ye_beq);
    saved.store(&cold, &eos, direct, rho, press, temp, Ye, thermal, use_press,
                Ye_beq);
    if (!load())
      CCTK_ERROR("Matching uniform atmosphere was not reused");
    CheckAtmoCopy(state, direct);

    tol = 0.02;
    if (!load() || state.rho_cut != direct.rho_atmo * (1 + tol))
      CCTK_ERROR("Global atmosphere ignored a changed cutoff");
    CheckAtmoCopy(saved.atmo, direct);
    for (CCTK_REAL bad : {-1.0, std::numeric_limits<CCTK_REAL>::infinity(),
                           std::numeric_limits<CCTK_REAL>::quiet_NaN()}) {
      tol = bad;
      if (load())
        CCTK_ERROR("Global atmosphere accepted an invalid cutoff");
    }
    tol = 0.001;

    nr = 1.0;
    if (load())
      CCTK_ERROR("Global atmosphere ignored density grading");
    nr = 0.0;
    np = 1.0;
    if (load() != (!thermal || !use_press))
      CCTK_ERROR("Incorrect pressure-grading reuse decision");
    np = 0.0;
    nt = 1.0;
    if (load() != (!thermal || use_press))
      CCTK_ERROR("Incorrect temperature-grading reuse decision");
    nt = 0.0;

    for (CCTK_REAL *value : {&rho, &Ye}) {
      *value *= 2;
      if (load())
        CCTK_ERROR("Global atmosphere accepted changed density or composition");
      *value /= 2;
    }
    press *= 2;
    if (load() != (!thermal || !use_press))
      CCTK_ERROR("Incorrect pressure-input reuse decision");
    press /= 2;
    temp *= 2;
    if (load() != (!thermal || use_press))
      CCTK_ERROR("Incorrect temperature-input reuse decision");
    temp /= 2;
    thermal = !thermal;
    if (load())
      CCTK_ERROR("Global atmosphere accepted a changed thermal mode");
    thermal = !thermal;
    use_press = !use_press;
    if (load() != !thermal)
      CCTK_ERROR("Incorrect pressure-mode reuse decision");
    use_press = !use_press;
    Ye_beq = true;
    if (load())
      CCTK_ERROR("Global atmosphere accepted a changed Ye mode");
    Ye_beq = false;

    eos_ptr = &other_eos;
    if (load())
      CCTK_ERROR("Global atmosphere accepted a different evolution EOS");
    eos_ptr = &eos;
    cold_ptr = &other_cold;
    if (load() != thermal)
      CCTK_ERROR("Incorrect cold-EOS reuse decision");
    cold_ptr = &cold;

    // Reinitialization must replace both the state and its input snapshot.
    rho *= 2;
    const auto updated =
        make_atmo(&cold, &eos, 0.0, rho, press, temp, Ye, 10.0,
                  nr, np, nt, tol, thermal, use_press, Ye_beq);
    saved.store(&cold, &eos, updated, rho, press, temp, Ye, thermal, use_press,
                Ye_beq);
    if (!load())
      CCTK_ERROR("Global atmosphere reinitialization failed");
    CheckAtmoCopy(state, updated);
  }
}

template <typename EOSIDType, typename EOSType>
void TestAtmoDevice(const EOSIDType *eos_1p, const EOSType *eos_3p) {
  DECLARE_CCTK_PARAMETERS;

  atmosphere state{};
  const bool graded = n_rho_atmo != 0.0 ||
      (thermal_eos_atmo &&
       (use_press_atmo ? n_press_atmo : n_temp_atmo) != 0.0);
  const bool use_global = get_global_atmo(eos_1p, eos_3p, state);
  if (use_global == graded)
    CCTK_ERROR("Incorrect global atmosphere selection for active parameters");
  if (use_global)
    CheckAtmoCopy(state, global_atmo.atmo);

  // Compare the startup state with device construction using the active EOS.
  // Zero exponents make the inner state valid at every supplied radius.
  const atmosphere expected = global_atmo.atmo;
  // Also exercise the C2P selection with active grading, across r_atmo.
  amrex::GpuArray<atmosphere, 8> atmo_ref{};
  for (int i = 0; i < 8; ++i) {
    const CCTK_REAL radial_distance = 0.5 * CCTK_REAL(i) * r_atmo;
    atmo_ref[i] = make_atmo(
        eos_1p, eos_3p, radial_distance, rho_abs_min, p_atmo, t_atmo,
        Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo, n_temp_atmo,
        atmo_tol, thermal_eos_atmo, use_press_atmo, Ye_atmo_beq);
  }
  // Neighboring centers straddle r_atmo for a face with spacing r_atmo/2.
  amrex::GpuArray<atmosphere, 2> face_ref{};
  for (int f = 0; f < 2; ++f) {
    face_ref[f] = make_atmo(
        eos_1p, eos_3p, (0.75 + 0.5 * f) * r_atmo,
        rho_abs_min, p_atmo, t_atmo, Ye_atmo, r_atmo,
        n_rho_atmo, n_press_atmo, n_temp_atmo, atmo_tol,
        thermal_eos_atmo, use_press_atmo, Ye_atmo_beq);
  }
  amrex::Gpu::DeviceScalar<unsigned int> failures(0);
  auto *failed = failures.dataPtr();
  amrex::ParallelFor(22, [=] AMREX_GPU_DEVICE(int i) {
    atmosphere actual{}, reference{};
    if (i < 8) {
      actual = make_atmo(
          eos_1p, eos_3p, CCTK_REAL(i), rho_abs_min, p_atmo, t_atmo, Ye_atmo,
          r_atmo, 0.0, 0.0, 0.0, atmo_tol, thermal_eos_atmo,
          use_press_atmo, Ye_atmo_beq);
      reference = expected;
    } else if (i < 16) {
      const CCTK_REAL radial_distance = 0.5 * CCTK_REAL(i - 8) * r_atmo;
      actual = use_global ? state : make_atmo(
          eos_1p, eos_3p, radial_distance, rho_abs_min, p_atmo, t_atmo,
          Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo, n_temp_atmo,
          atmo_tol, thermal_eos_atmo, use_press_atmo, Ye_atmo_beq);
      reference = atmo_ref[i - 8];
    } else {
      // Check both neighboring positions for x-, y- and z-directed faces.
      const int dir_i = (i - 16) / 2;
      const int f = (i - 16) % 2;
      CCTK_REAL r2_atm = 0.0;
      for (int ii = 0; ii < 3; ++ii) {
        const CCTK_REAL x = ii == dir_i ? r_atmo : 0.0;
        const CCTK_REAL dx = (ii == dir_i) * 0.5 * (0.5 * r_atmo);
        r2_atm += (x + (2 * f - 1) * dx) * (x + (2 * f - 1) * dx);
      }
      actual = use_global ? state : make_atmo(
          eos_1p, eos_3p, sqrt(r2_atm), rho_abs_min, p_atmo, t_atmo,
          Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo, n_temp_atmo,
          atmo_tol, thermal_eos_atmo, use_press_atmo, Ye_atmo_beq);
      reference = face_ref[f];
    }
    const CCTK_REAL a[] = {actual.rho_atmo, actual.eps_atmo, actual.ye_atmo,
                          actual.press_atmo, actual.temp_atmo,
                          actual.entropy_atmo, actual.rho_cut};
    const CCTK_REAL b[] = {reference.rho_atmo, reference.eps_atmo, reference.ye_atmo,
                          reference.press_atmo, reference.temp_atmo,
                          reference.entropy_atmo, reference.rho_cut};
    for (int n = 0; n < 7; ++n) {
      const CCTK_REAL scale = fmax(fabs(b[n]), n == 5 ? 1.0 : 1.0e-12);
      if (!std::isfinite(a[n]) || !std::isfinite(b[n]) ||
          fabs(a[n] - b[n]) > 1.0e-7 * scale)
        amrex::HostDevice::Atomic::Add(failed, 1U);
    }
  });
  amrex::Gpu::streamSynchronize();
  if (failures.dataValue())
    CCTK_ERROR("Global atmosphere differs from device construction");
}
} // namespace

extern "C" void AsterX_TestAtmo(CCTK_ARGUMENTS) {
  TestAtmoGlobal();
  if (global_atmo.valid) {
    const auto test = [](const auto *eos_3p) {
      if (global_eos_1p_pwpoly)
        TestAtmoDevice(global_eos_1p_pwpoly, eos_3p);
      else
        TestAtmoDevice(global_eos_1p_poly, eos_3p);
    };
    if (global_eos_3p_ig)
      test(global_eos_3p_ig);
    else if (global_eos_3p_tab3d)
      test(global_eos_3p_tab3d);
  }
  CCTK_INFO("Global atmosphere tests passed");
}

} // namespace AsterX
