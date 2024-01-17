#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "c2p_utils.hxx"
#include "c2p.hxx"
#include "c2p_1DPalenzuela.hxx"
#include "c2p_2DNoble.hxx"
#include <eos.hxx>
#include <eos_idealgas.hxx>

namespace Con2PrimFactory {

using namespace Arith;
using namespace EOSX;

extern "C" void Con2PrimFactory_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  // Initializing the IG EOS
  eos::range rgeps(0, 100), rgrho(1e-13, 20), rgye(0.5, 0.5);
  const eos_idealgas eos_th(2.0, 938.985, rgeps, rgrho, rgye);

  // Setting up atmosphere
  atmosphere atmo(1e-10, 1e-8, 0.5, 1e-8, 0.001);

  // Metric
  const smat<CCTK_REAL, 3> g{1.0, 0.0, 0.0,
                             1.0, 0.0, 1.0}; // xx, xy, xz, yy, yz, zz

  // Con2Prim objects
  c2p_2DNoble c2p_Noble(eos_th, atmo, 100, 1e-8, 1e8, 1, 1, true);
  c2p_1DPalenzuela c2p_Pal(eos_th, atmo, 100, 1e-8, 1e8, 1, 1, true);

  // Construct error report object:
  c2p_report rep_Noble;
  c2p_report rep_Pal;

  prim_vars pv;
  // rho(p.I), eps(p.I), dummy_Ye, press(p.I),v_up, wlor, Bup
  prim_vars pv_seeds{0.125, 0.8, 0.5, 0.1, {0, 0, 0}, 1, {0.5, -1.0, 0}};

  // cons_vars cv{dens(p.I), {momx(p.I), momy(p.I), momz(p.I)}, tau(p.I),
  //  dummy_dYe, {dBx(p.I), dBy(p.I), dBz(p.I)}};

  cons_vars cv;
  cv.from_prim(pv_seeds, g);

  // Testing C2P Noble
  CCTK_VINFO("Testing C2P Noble...");
  c2p_Noble.solve(eos_th, pv, pv_seeds, cv, g, rep_Noble);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.vel(0), pv.vel(0), pv_seeds.vel(1),
         pv.vel(1), pv_seeds.vel(2), pv.vel(2), pv_seeds.Bvec(0), pv.Bvec(0),
         pv_seeds.Bvec(1), pv.Bvec(1), pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "dYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n",
         cv.dens, cv.tau, cv.mom(0), cv.mom(1), cv.mom(2), cv.dYe, cv.dBvec(0),
         cv.dBvec(1), cv.dBvec(2));
  /*
    assert(pv.rho == pv_seeds.rho);
    assert(pv.eps == pv_seeds.eps);
    assert(pv.press == pv_seeds.press);
    assert(pv.vel == pv_seeds.vel);
    assert(pv.Bvec == pv_seeds.Bvec);
  */

  rep_Noble.debug_message();

  // Testing C2P Palenzuela
  CCTK_VINFO("Testing C2P Palenzuela...");
  c2p_Pal.solve(eos_th, pv, pv_seeds, cv, g, rep_Pal);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.vel(0), pv.vel(0), pv_seeds.vel(1),
         pv.vel(1), pv_seeds.vel(2), pv.vel(2), pv_seeds.Bvec(0), pv.Bvec(0),
         pv_seeds.Bvec(1), pv.Bvec(1), pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "dYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n",
         cv.dens, cv.tau, cv.mom(0), cv.mom(1), cv.mom(2), cv.dYe, cv.dBvec(0),
         cv.dBvec(1), cv.dBvec(2));
  /*
    assert(pv.rho == pv_seeds.rho);
    assert(pv.eps == pv_seeds.eps);
    assert(pv.press == pv_seeds.press);
    assert(pv.vel == pv_seeds.vel);
    assert(pv.Bvec == pv_seeds.Bvec);
  */
  rep_Pal.debug_message();
}

} // namespace Con2PrimFactory
