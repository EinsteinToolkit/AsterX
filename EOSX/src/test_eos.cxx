#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "setup_eos.hxx"

namespace EOSX {

// using namespace Arith;
// using namespace Con2PrimFactory;
// using namespace AsterUtils;

extern "C" void EOS_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;


  // Get local eos objects
  auto eos_3p_hyb = global_eos_3p_hyb_pwpoly;
  auto eos_1p_pwpoly = global_eos_1p_pwpoly;
  
  //pick test value(s) of rho
  const CCTK_REAL rho_test = 3.15519e-04;

  //compute pressure to see if it matches what we expect
  const CCTK_REAL press_test_hyb =
      eos_3p_hyb->p_cold(rho_test);

  // print cold pressure via hybrid EOS
  // Testing Peicewise polytropic EOS
  CCTK_VINFO("Testing PWpoly EOS...");

  printf("rho_test: %f \n"
         "press: %f \n",
         rho_test, press_test_hyb);
  /*
    assert(pv.rho == pv_seeds.rho);
  */


  ////////////////////////////////////////////
  //check consitency of segment functions  //
  //////////////////////////////////////////

  // compute enthalpy from rho
  const CCTK_REAL gm1_test = 
      eos_1p_pwpoly->segment_for_rho(rho_test).gm1_from_rho(rho_test);

  // compute pressure from enthalpy using rho & enthalpy to find segment
  const CCTK_REAL press_test_rho_gm1 =
      eos_1p_pwpoly->segment_for_rho(rho_test).p_from_gm1(gm1_test);

  const CCTK_REAL press_test_gm1_gm1 =
      eos_1p_pwpoly->segment_for_gm1(gm1_test).p_from_gm1(gm1_test);

  // recompute density from enthalpy using all three meathods to find segment
  const CCTK_REAL recomp_rho_test_rho_gm1 =
      eos_1p_pwpoly->segment_for_rho(rho_test).rho_from_gm1(gm1_test);

  const CCTK_REAL recomp_rho_test_gm1_gm1 =
      eos_1p_pwpoly->segment_for_gm1(gm1_test).rho_from_gm1(gm1_test);

  const CCTK_REAL recomp_rho_test_press_gm1 =
      eos_1p_pwpoly->segment_for_p(press_test_rho_gm1).rho_from_gm1(gm1_test);

  // print pressures and recomputed densities using different segment finders
  // Testing Peicewise polytropic EOS
  CCTK_VINFO("Testing PWpoly EOS...");

  printf("Pressures using different segments: %f, %f \n"
         "Densities using different segments: %f, %f, %f, %f \n",
         press_test_rho_gm1, press_test_gm1_gm1, rho_test,
         recomp_rho_test_rho_gm1, recomp_rho_test_gm1_gm1,
         recomp_rho_test_press_gm1);
  /*
    assert(pv.rho == pv_seeds.rho);
  */

  //other ideas
  //recompute rho from second pressure?
  //recompute enthalpy from pressures?
  // auto seg_rho_test =
  //     eos_1p_pwpoly->segment_for_rho(rho_test);
  // auto seg_press_test =
  //     eos_1p_pwpoly->segment_for_press(press_test_hyb);



//ALF2 Params
// #atmosphere for Ideal gas EOS
// Con2PrimFactory::thermal_eos_atmo = "no"
// Con2PrimFactory::rho_abs_min = 1e-11
// Con2PrimFactory::atmo_tol = 1e-3
// Con2PrimFactory::t_atmo = 9.9770006382257168e-03 #0.05
// Con2PrimFactory::Ye_atmo = 0.5
// #Equation of State
// #---------------------------------
// EOSX::evolution_eos = "Hybrid"
// #EOSX::gl_gamma = 2.0
// #EOSX::poly_gamma = 2.0
// #EOSX::poly_k = 123.6
// EOSX::rho_max = 1e8
// EOSX::eps_max = 1e8
// EOSX::eps_min = 1.236e-09

// #ALF2 parameters
// #PW Polytropic params
// EOSX::pwpoly_nsegm = 4
// EOSX::pwpoly_rho_p0 = 865.203 #need to double check
// EOSX::pwpoly_segm_bound = [3.15519e-04, 8.11939e-04, 1.62003e-03]
// EOSX::pwpoly_segm_gamma = [1.35692, 4.07, 2.411, 1.89]

// #Hybrid params
// EOSX::gamma_th = 1.6666666666666667 #1.8


}

} // namespace EOSX
