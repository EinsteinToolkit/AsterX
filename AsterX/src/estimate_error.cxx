#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <cmath>

namespace AsterX {

/* calculate max of d(var)/(dx) * hx in all dirs */
template <typename T>
CCTK_DEVICE CCTK_HOST inline T calc_grad_1st(const Loop::GF3D2<const T> &gf,
                                             const Loop::PointDesc &p) {
  using std::max, std::fabs;

  CCTK_REAL err{0}, errp{0}, errm{0};

  for (int d = 0; d < Loop::dim; ++d) {
    auto varm = gf(p.I - p.DI[d]);
    auto var0 = gf(p.I);
    auto varp = gf(p.I + p.DI[d]);
    errp += (varp - var0) * (varp - var0);
    errm += (var0 - varm) * (var0 - varm);
  }

  err = max({errp, errm});

  return sqrt(err);
}

/* calculate max of d(var)/(dx) * hx in all dirs */
template <typename T>
CCTK_DEVICE CCTK_HOST inline T calc_deriv_1st(const Loop::GF3D2<const T> &gf,
                                              const Loop::PointDesc &p) {
  using std::max, std::fabs;

  CCTK_REAL err{0};

  for (int d = 0; d < Loop::dim; ++d) {
    auto varm = gf(p.I - p.DI[d]);
    auto var0 = gf(p.I);
    auto varp = gf(p.I + p.DI[d]);
    err = max({err, fabs(var0 - varm), fabs(varp - var0)});
  }

  return err;
}

/* calculate max of d^2(var)/(dx^2) * hx^2 in all dirs */
template <typename T>
CCTK_DEVICE CCTK_HOST inline T calc_deriv_2nd(const Loop::GF3D2<const T> &gf,
                                              const Loop::PointDesc &p) {
  using std::max, std::fabs;

  CCTK_REAL err{0};

  for (int d = 0; d < Loop::dim; ++d) {
    auto varm = gf(p.I - p.DI[d]);
    auto var0 = gf(p.I);
    auto varp = gf(p.I + p.DI[d]);
    err = max({err, fabs(varp + varm - 2 * var0)});
  }

  return err;
}

/* calculate max of second derivative error norm in all dirs
 * based on eq. 50 of Mignone+2011 (https://arxiv.org/abs/1110.0740) */
template <typename T>
CCTK_DEVICE CCTK_HOST inline T
calc_deriv_2nd_norm(const Loop::GF3D2<const T> &gf, const Loop::PointDesc &p,
                    CCTK_REAL epsilon_err) {
  using std::max, std::fabs, std::sqrt;

  CCTK_REAL err{0};

  for (int d = 0; d < Loop::dim; ++d) {
    auto varm = gf(p.I - p.DI[d]);
    auto var0 = gf(p.I);
    auto varp = gf(p.I + p.DI[d]);

    auto diffp0 = varp - var0;
    auto diff0m = var0 - varm;

    auto varsigma = fabs(varp) + 2 * fabs(var0) + fabs(varm);

    auto diffpm2 = fabs(diffp0 - diff0m) * fabs(diffp0 - diff0m);
    auto diffpm_denom = (fabs(diffp0) + fabs(diff0m) + epsilon_err * varsigma) *
                        (fabs(diffp0) + fabs(diff0m) + epsilon_err * varsigma);

    auto chi = sqrt(diffpm2 / diffpm_denom);
    err = max({err, chi});
  }

  return err;
}

extern "C" void AsterX_Error_fd_rho(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Error_fd_rho;
  using std::max;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        regrid_error(p.I) = max(0.0, calc_deriv_1st(rho, p));
      });
}

extern "C" void AsterX_Error_sd_rho(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Error_sd_rho;

  using std::max;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        regrid_error(p.I) = max(0.0, calc_deriv_2nd(rho, p));
      });
}

extern "C" void AsterX_Error_sdn_rho(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Error_sdn_rho;
  DECLARE_CCTK_PARAMETERS;

  using std::max;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        regrid_error(p.I) = max(0.0, calc_deriv_2nd_norm(rho, p, epsilon_err));
      });
}

extern "C" void AsterX_Error_fg_rho(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Error_fg_rho;
  DECLARE_CCTK_PARAMETERS;

  using std::max;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        regrid_error(p.I) = max(0.0, calc_grad_1st(rho, p));
      });
}

} // namespace AsterX
