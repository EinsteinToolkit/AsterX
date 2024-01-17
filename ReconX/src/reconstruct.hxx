#ifndef RECONX_RECONSTRUCT_HXX
#define RECONX_RECONSTRUCT_HXX

#include <loop_device.hxx>

#include <cctk.h>
// #include <cctk_Arguments.h>
// #include <cctk_Parameters.h>

#include "monocentral.hxx"
#include "minmod.hxx"
#include "ppm.hxx"
#include "eppm.hxx"
#include "wenoz.hxx"
#include "mp5.hxx"

#include <array>

namespace ReconX {

using std::array;
using namespace Arith;

// enum class for different reconstruction routines

enum class reconstruction_t {
  Godunov,
  minmod,
  monocentral,
  ppm,
  eppm,
  wenoz,
  mp5
};

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE array<CCTK_REAL, 2>
reconstruct(const GF3D2<const CCTK_REAL> &gf_var, const PointDesc &p,
            const reconstruction_t &reconstruction, const int &dir,
            const bool &gf_is_rho, const bool &gf_is_press,
            const GF3D2<const CCTK_REAL> &gf_press,
            const GF3D2<const CCTK_REAL> &gf_vel_dir,
            const reconstruct_params_t &reconstruct_params) {
  // Neighbouring "plus" and "minus" cell indices
  const auto Immm = p.I - 3 * p.DI[dir];
  const auto Imm = p.I - 2 * p.DI[dir];
  const auto Im = p.I - p.DI[dir];
  const auto Ip = p.I;
  const auto Ipp = p.I + p.DI[dir];
  const auto Ippp = p.I + 2 * p.DI[dir];

  switch (reconstruction) {

  case reconstruction_t::Godunov: {
    return {gf_var(Im), gf_var(Ip)};
  }

  case reconstruction_t::minmod: {
    return minmod_reconstruct(gf_var(Imm), gf_var(Im), gf_var(Ip), gf_var(Ipp));
  }

  case reconstruction_t::monocentral: {
    return monocentral_reconstruct(gf_var(Imm), gf_var(Im), gf_var(Ip),
                                   gf_var(Ipp));
  }

  case reconstruction_t::ppm: {
    return ppm_reconstruct(
        gf_var(Immm), gf_var(Imm), gf_var(Im), gf_var(Ip), gf_var(Ipp),
        gf_var(Ippp), gf_press(Immm), gf_press(Imm), gf_press(Im), gf_press(Ip),
        gf_press(Ipp), gf_press(Ippp), gf_vel_dir(Imm), gf_vel_dir(Im),
        gf_vel_dir(Ip), gf_vel_dir(Ipp), gf_is_rho, reconstruct_params);
  }

  case reconstruction_t::wenoz: {
    return wenoz_reconstruct(gf_var(Immm), gf_var(Imm), gf_var(Im), gf_var(Ip),
                             gf_var(Ipp), gf_var(Ippp),
                             reconstruct_params.weno_eps);
  }

  case reconstruction_t::mp5: {
    return mp5_reconstruct(gf_var(Immm), gf_var(Imm), gf_var(Im), gf_var(Ip),
                           gf_var(Ipp), gf_var(Ippp),
                           reconstruct_params.mp5_alpha);
  }

  case reconstruction_t::eppm: {
    const array<const vect<int, dim>, 5> cells_Im = {Immm, Imm, Im, Ip, Ipp};
    const array<const vect<int, dim>, 5> cells_Ip = {Imm, Im, Ip, Ipp, Ippp};

    const array<CCTK_REAL, 2> rc_Im =
        eppm(gf_var, cells_Im, gf_is_press, gf_press, gf_vel_dir,
             reconstruct_params);
    const array<CCTK_REAL, 2> rc_Ip =
        eppm(gf_var, cells_Ip, gf_is_press, gf_press, gf_vel_dir,
             reconstruct_params);

    return array<CCTK_REAL, 2>{rc_Im[1], rc_Ip[0]};
  }

  default:
    assert(0);
  }
}

} // namespace ReconX

#endif // RECONX_RECONSTRUCT_HXX
