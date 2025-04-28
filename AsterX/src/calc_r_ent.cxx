#include <loop_device.hxx>

#include <derivs.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

#include "aster_utils.hxx"
#include "eos.hxx"

namespace AsterX {
using namespace Loop;
using namespace Arith;
using namespace AsterUtils;

enum class eos_t { IdealGas, Hybrid, Tabulated };

extern "C" void AsterX_CalcPhysEntropy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CalcPhysEntropy;
  DECLARE_CCTK_PARAMETERS;

  // First calculate "physical" entropy from "evolved" entropy
  // defining EOS objects
  eos_t eostype;
  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eostype = eos_t::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eostype = eos_t::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated")) {
    eostype = eos_t::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        switch (eostype) {
        case eos_t::IdealGas: {
          phys_ent(p.I) = std::log(entropy(p.I));
          break;
        }
        case eos_t::Hybrid: {
          printf("Hybrid EOS is not yet supported");
          assert(0);
          break;
        }
        case eos_t::Tabulated: {
          printf("Tabulated EOS is not yet supported");
          assert(0);
          break;
        }
        default:
          assert(0);
        }
      });
}

extern "C" void AsterX_CalcEntropyResidual(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CalcEntropyResidual;
  DECLARE_CCTK_PARAMETERS;

  // Compute d_i (phys entropy)
  // Derivs boilerplate
  Arith::vect<int, dim> imin, imax;
  const std::array<int, dim> nghostzones = {
      cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]};
  GridDescBase(cctkGH).box_int<1, 1, 1>(nghostzones, imin, imax);
  const GF3D5layout layout5(imin, imax);
  const Arith::vect<CCTK_REAL, dim> dx(std::array<CCTK_REAL, dim>{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  });

  constexpr int nvars = 4;
  int ivar = 0;
  GF3D5vector<CCTK_REAL> vars(layout5, nvars);
  const auto make_gf = [&]() { return GF3D5<CCTK_REAL>(vars(ivar++)); };
  const auto make_vec = [&](const auto &f) {
    return Arith::vec<std::result_of_t<decltype(f)()>, dim>(
        [&](int) { return f(); });
  };
  const auto make_vec_gf = [&]() { return make_vec(make_gf); };

  const GF3D5<CCTK_REAL> t5_s(make_gf());
  const Arith::vec<GF3D5<CCTK_REAL>, dim> t5_ds(make_vec_gf());

  Derivs::calc_derivs<1, 1, 1>(t5_s, t5_ds, layout5, grid, phys_ent, dx,
                               efl_deriv_order);

  // Prep spacetime GFs
  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Calculate spatial contraction
        const GF3D5index index5(layout5, p.I);
        const vec<CCTK_REAL, 3> di_s = t5_ds(index5);

        const CCTK_REAL alp_avg = calc_avg_v2c(alp, p);
        const vec<CCTK_REAL, 3> betas_avg(
            [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p); });
        const vec<CCTK_REAL, 3> vels{velx(p.I), vely(p.I), velx(p.I)};

        const CCTK_REAL v_dis =
            calc_contraction(alp_avg * vels - betas_avg, di_s);

        // Calculate d_t s
        const CCTK_REAL i2dt = 1 / (2 * cctk_delta_time);
        const CCTK_REAL dts =
            i2dt * (3 * phys_ent(p.I) - 4 * ent_m1(p.I) + ent_m2(p.I));

        // Calculate R
        r_ent(p.I) = std::abs(dts + v_dis);

        // Cycle timelevels
        ent_m2(p.I) = ent_m1(p.I);
        ent_m1(p.I) = phys_ent(p.I);
      });
}

} // namespace AsterX
