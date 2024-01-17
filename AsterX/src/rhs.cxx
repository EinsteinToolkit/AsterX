#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "utils.hxx"

// #ifdef AMREX_USE_GPU
// #include <AMReX_GpuDevice.H>
// #endif

#include <array>

namespace AsterX {
using namespace Loop;
using namespace Arith;

enum class vector_potential_gauge_t { algebraic, generalized_lorentz };

extern "C" void AsterX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_RHS;
  DECLARE_CCTK_PARAMETERS;

  vector_potential_gauge_t gauge;
  if (CCTK_EQUALS(vector_potential_gauge, "algebraic"))
    gauge = vector_potential_gauge_t::algebraic;
  else if (CCTK_EQUALS(vector_potential_gauge, "generalized Lorentz"))
    gauge = vector_potential_gauge_t::generalized_lorentz;
  else
    CCTK_ERROR("Unknown value for parameter \"vector_potential_gauge\"");

  const vec<CCTK_REAL, dim> idx{1 / CCTK_DELTA_SPACE(0),
                                1 / CCTK_DELTA_SPACE(1),
                                1 / CCTK_DELTA_SPACE(2)};

  const vec<GF3D2<const CCTK_REAL>, dim> gf_fdens{fxdens, fydens, fzdens};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_fmomx{fxmomx, fymomx, fzmomx};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_fmomy{fxmomy, fymomy, fzmomy};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_fmomz{fxmomz, fymomz, fzmomz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_ftau{fxtau, fytau, fztau};
  const vec<vec<GF3D2<const CCTK_REAL>, dim>, dim> gf_fBs{
      {fxBx, fyBx, fzBx}, {fxBy, fyBy, fzBy}, {fxBz, fyBz, fzBz}};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_F{Fx, Fy, Fz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_Fbeta{Fbetax, Fbetay, Fbetaz};

  const auto calcupdate_hydro =
      [=] CCTK_DEVICE(const vec<GF3D2<const CCTK_REAL>, dim> &gf_fluxes,
                      const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        vec<CCTK_REAL, 3> dfluxes([&](int i) ARITH_INLINE {
          return gf_fluxes(i)(p.I + p.DI[i]) - gf_fluxes(i)(p.I);
        });
        return -calc_contraction(idx, dfluxes);
      };

  const auto calcupdate_Avec =
      [=] CCTK_DEVICE(const PointDesc &p, int i) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int j = (i == 0) ? 1 : ((i == 1) ? 2 : 0);
        const int k = (i == 0) ? 2 : ((i == 1) ? 0 : 1);

        const CCTK_REAL E =
            0.25 * ((gf_fBs(j)(k)(p.I) + gf_fBs(j)(k)(p.I - p.DI[j])) -
                    (gf_fBs(k)(j)(p.I) + gf_fBs(k)(j)(p.I - p.DI[k])));

        switch (gauge) {
        case vector_potential_gauge_t::algebraic: {
          return -E;
          break;
        }
        case vector_potential_gauge_t::generalized_lorentz: {
          return -E - calc_fd2_v2e(G, p, i);
          break;
        }
        default:
          assert(0);
        }
      };

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        densrhs(p.I) += calcupdate_hydro(gf_fdens, p);
        momxrhs(p.I) += calcupdate_hydro(gf_fmomx, p);
        momyrhs(p.I) += calcupdate_hydro(gf_fmomy, p);
        momzrhs(p.I) += calcupdate_hydro(gf_fmomz, p);
        taurhs(p.I) += calcupdate_hydro(gf_ftau, p);
        if (isnan(densrhs(p.I))) {
          printf("calcupdate = %f, ", calcupdate_hydro(gf_fdens, p));
          printf("densrhs = %f, gf_fdens = %f, %f, %f, %f, %f, %f \n",
              densrhs(p.I),
              gf_fdens(0)(p.I), gf_fdens(1)(p.I), gf_fdens(2)(p.I),
              gf_fdens(0)(p.I + p.DI[0]), gf_fdens(1)(p.I + p.DI[1]), gf_fdens(2)(p.I + p.DI[2]));
        }
        assert(!isnan(densrhs(p.I)));
      });

  grid.loop_int_device<1, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      Avec_x_rhs(p.I) = calcupdate_Avec(p, 0);
                                    });

  grid.loop_int_device<0, 1, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      Avec_y_rhs(p.I) = calcupdate_Avec(p, 1);
                                    });

  grid.loop_int_device<0, 0, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      Avec_z_rhs(p.I) = calcupdate_Avec(p, 2);
                                    });

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        switch (gauge) {
        case vector_potential_gauge_t::algebraic: {
          Psi_rhs(p.I) = 0.0;
          break;
        }

        case vector_potential_gauge_t::generalized_lorentz: {
          CCTK_REAL dF = 0.0;
          for (int i = 0; i < dim; i++) {
            /* diFi on vertices (should be v2v but c2c works too) */
            dF += calc_fd2_c2c(gf_F(i), p, i) -
                  (gf_beta(i)(p.I) < 0
                       ? calc_fd2_v2v_oneside(gf_Fbeta(i), p, i, -1)
                       : calc_fd2_v2v_oneside(gf_Fbeta(i), p, i, 1));
          }
          Psi_rhs(p.I) = -dF - lorenz_damp_fac * alp(p.I) * Psi(p.I);
          break;
        }

        default:
          assert(0);
        }
      });
}

} // namespace AsterX
