#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <seeds_utils.hxx>

#include <eos.hxx>
#include <eos_idealgas.hxx>

namespace AsterSeeds {
using namespace std;
using namespace Loop;
using namespace EOSX;

extern "C" void Tests1D_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Tests1D_Initialize;
  DECLARE_CCTK_PARAMETERS;

  // For all the tests, the initial data EOS is ideal gas
  // Constructing the IG EOS object
  eos::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
      rgye(ye_min, ye_max);

  const eos_idealgas eos_th(gl_gamma, particle_mass, rgeps, rgrho, rgye);
  const CCTK_REAL dummy_ye = 0.5;

  if (CCTK_EQUALS(test_case, "equilibrium")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0;
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
          press(p.I) = 1.0;
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

  } else if (CCTK_EQUALS(test_case, "sound wave")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0 + amplitude * sin(M_PI * p.x);
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
          press(p.I) = 1.0; // should add kinetic energy here
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });
  
  } else if (CCTK_EQUALS(test_case, "Alfven wave")) {
    const CCTK_REAL A0 = 1.0;
    const CCTK_REAL va = 0.5;
    const CCTK_REAL k = 2*M_PI;

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0;
          vely(p.I) = -va * A0 * cos(k*p.x);
          velz(p.I) = -va * A0 * sin(k*p.x);
          press(p.I) = 0.5; // should add kinetic energy here
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = p.z*cos(k*p.x) - p.y*sin(k*p.x); });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = -p.z/2.0; });
            //CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = p.x*sin(k*p.x) - p.z; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = p.y/2.0; });
            //CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = p.y - p.x*cos(k*p.x); });

  } else if (CCTK_EQUALS(test_case, "shock tube")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = 2.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 2.0;
          } else {
            rho(p.I) = 1.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 1.0;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

  } else if (CCTK_EQUALS(test_case, "Balsara1")) {
   
    CCTK_REAL rhol = 1.0;
    CCTK_REAL vxl = 0.0;
    CCTK_REAL vyl = 0.0;
    CCTK_REAL vzl = 0.0;
    CCTK_REAL pressl = 1.0;
    CCTK_REAL Bxl = 0.5;
    CCTK_REAL Byl = 1.0;
    CCTK_REAL Bzl = 0.0;

    CCTK_REAL rhor = 0.125;
    CCTK_REAL vxr = 0.0;
    CCTK_REAL vyr = 0.0;
    CCTK_REAL vzr = 0.0;
    CCTK_REAL pressr = 0.1;
    CCTK_REAL Bxr = 0.5;
    CCTK_REAL Byr = -1.0;
    CCTK_REAL Bzr = 0.0;    
	 
    if (CCTK_EQUALS(shock_dir, "x")) {
   
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = rhol;
            velx(p.I) = vxl;
            vely(p.I) = vyl;
            velz(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velx(p.I) = vxr;
            vely(p.I) = vyr;
            velz(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.x <= 0.0) { 
	       Avec_x(p.I) = Byl * (p.z) - Bzl * (p.y);  
	    } else { 
	       Avec_x(p.I) = Byr * (p.z) - Bzr * (p.y); } });


    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_z(p.I) = Bxr * (p.y);  });
    
    } else if (CCTK_EQUALS(shock_dir, "y")) {

    //x-->y, y-->z, z-->x
    Byl = 0.5;
    Bzl = 1.0;
    Bxl = 0.0;
    Byr = 0.5;
    Bzr = -1.0;
    Bxr = 0.0;

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.y <= 0.0) {
            rho(p.I) = rhol;
            vely(p.I) = vxl;
            velz(p.I) = vyl;
            velx(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            vely(p.I) = vxr;
            velz(p.I) = vyr;
            velx(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.y <= 0.0) {
               Avec_y(p.I) = Bzl * (p.x) - Bxl * (p.z);
            } else {
               Avec_y(p.I) = Bzr * (p.x) - Bxr * (p.z); } });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_x(p.I) = Byr * (p.z);  });

    } else if (CCTK_EQUALS(shock_dir, "z")) {
    
    //x-->z, y-->x, z-->y
    Bzl = 0.5;
    Bxl = 1.0;
    Byl = 0.0;
    Bzr = 0.5;
    Bxr = -1.0;
    Byr = 0.0;

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.z <= 0.0) {
            rho(p.I) = rhol;
            velz(p.I) = vxl;
            velx(p.I) = vyl;
            vely(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velz(p.I) = vxr;
            velx(p.I) = vyr;
            vely(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.z <= 0.0) {
               Avec_z(p.I) = Bxl * (p.y) - Byl * (p.x);
            } else {
               Avec_z(p.I) = Bxr * (p.y) - Byr * (p.x); } });


    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_y(p.I) = Bzr * (p.x);  });

    } else { CCTK_ERROR("Shock direction case not defined"); }
    
  } else if (CCTK_EQUALS(test_case, "Balsara2")) {
   
    CCTK_REAL rhol = 1.0;
    CCTK_REAL vxl = 0.0;
    CCTK_REAL vyl = 0.0;
    CCTK_REAL vzl = 0.0;
    CCTK_REAL pressl = 30.0;
    CCTK_REAL Bxl = 5.0;
    CCTK_REAL Byl = 6.0;
    CCTK_REAL Bzl = 6.0;

    CCTK_REAL rhor = 1.0;
    CCTK_REAL vxr = 0.0;
    CCTK_REAL vyr = 0.0;
    CCTK_REAL vzr = 0.0;
    CCTK_REAL pressr = 1.0;
    CCTK_REAL Bxr = 5.0;
    CCTK_REAL Byr = 0.7;
    CCTK_REAL Bzr = 0.7;    
	 
    if (CCTK_EQUALS(shock_dir, "x")) {
   
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = rhol;
            velx(p.I) = vxl;
            vely(p.I) = vyl;
            velz(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velx(p.I) = vxr;
            vely(p.I) = vyr;
            velz(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.x <= 0.0) { 
	       Avec_x(p.I) = Byl * (p.z) - Bzl * (p.y);  
	    } else { 
	       Avec_x(p.I) = Byr * (p.z) - Bzr * (p.y); } });


    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_z(p.I) = Bxr * (p.y);  });
    
    } else if (CCTK_EQUALS(shock_dir, "y")) {

    //x-->y, y-->z, z-->x
    Byl = 5.0;
    Bzl = 6.0;
    Bxl = 6.0;
    Byr = 5.0;
    Bzr = 0.7;
    Bxr = 0.7;

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.y <= 0.0) {
            rho(p.I) = rhol;
            vely(p.I) = vxl;
            velz(p.I) = vyl;
            velx(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            vely(p.I) = vxr;
            velz(p.I) = vyr;
            velx(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.y <= 0.0) {
               Avec_y(p.I) = Bzl * (p.x) - Bxl * (p.z);
            } else {
               Avec_y(p.I) = Bzr * (p.x) - Bxr * (p.z); } });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_x(p.I) = Byr * (p.z);  });

    } else if (CCTK_EQUALS(shock_dir, "z")) {
    
    //x-->z, y-->x, z-->y
    Bzl = 5.0;
    Bxl = 6.0;
    Byl = 6.0;
    Bzr = 5.0;
    Bxr = 0.7;
    Byr = 0.7;

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.z <= 0.0) {
            rho(p.I) = rhol;
            velz(p.I) = vxl;
            velx(p.I) = vyl;
            vely(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velz(p.I) = vxr;
            velx(p.I) = vyr;
            vely(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.z <= 0.0) {
               Avec_z(p.I) = Bxl * (p.y) - Byl * (p.x);
            } else {
               Avec_z(p.I) = Bxr * (p.y) - Byr * (p.x); } });


    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_y(p.I) = Bzr * (p.x);  });

    } else { CCTK_ERROR("Shock direction case not defined"); }
    
  } else if (CCTK_EQUALS(test_case, "Balsara3")) {
   
    CCTK_REAL rhol = 1.0;
    CCTK_REAL vxl = 0.0;
    CCTK_REAL vyl = 0.0;
    CCTK_REAL vzl = 0.0;
    CCTK_REAL pressl = 1000.0;
    CCTK_REAL Bxl = 10.0;
    CCTK_REAL Byl = 7.0;
    CCTK_REAL Bzl = 7.0;

    CCTK_REAL rhor = 1.0;
    CCTK_REAL vxr = 0.0;
    CCTK_REAL vyr = 0.0;
    CCTK_REAL vzr = 0.0;
    CCTK_REAL pressr = 0.1;
    CCTK_REAL Bxr = 10.0;
    CCTK_REAL Byr = 0.7;
    CCTK_REAL Bzr = 0.7;    
	 
    if (CCTK_EQUALS(shock_dir, "x")) {
   
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = rhol;
            velx(p.I) = vxl;
            vely(p.I) = vyl;
            velz(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velx(p.I) = vxr;
            vely(p.I) = vyr;
            velz(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.x <= 0.0) { 
	       Avec_x(p.I) = Byl * (p.z) - Bzl * (p.y);  
	    } else { 
	       Avec_x(p.I) = Byr * (p.z) - Bzr * (p.y); } });


    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_z(p.I) = Bxr * (p.y);  });
    
    } else if (CCTK_EQUALS(shock_dir, "y")) {

    //x-->y, y-->z, z-->x
    Byl = 10.0;
    Bzl = 7.0;
    Bxl = 7.0;
    Byr = 10.0;
    Bzr = 0.7;
    Bxr = 0.7;

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.y <= 0.0) {
            rho(p.I) = rhol;
            vely(p.I) = vxl;
            velz(p.I) = vyl;
            velx(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            vely(p.I) = vxr;
            velz(p.I) = vyr;
            velx(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.y <= 0.0) {
               Avec_y(p.I) = Bzl * (p.x) - Bxl * (p.z);
            } else {
               Avec_y(p.I) = Bzr * (p.x) - Bxr * (p.z); } });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_x(p.I) = Byr * (p.z);  });

    } else if (CCTK_EQUALS(shock_dir, "z")) {
    
    //x-->z, y-->x, z-->y
    Bzl = 10.0;
    Bxl = 7.0;
    Byl = 7.0;
    Bzr = 10.0;
    Bxr = 0.7;
    Byr = 0.7;

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.z <= 0.0) {
            rho(p.I) = rhol;
            velz(p.I) = vxl;
            velx(p.I) = vyl;
            vely(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velz(p.I) = vxr;
            velx(p.I) = vyr;
            vely(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.z <= 0.0) {
               Avec_z(p.I) = Bxl * (p.y) - Byl * (p.x);
            } else {
               Avec_z(p.I) = Bxr * (p.y) - Byr * (p.x); } });


    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_y(p.I) = Bzr * (p.x);  });

    } else { CCTK_ERROR("Shock direction case not defined"); }
  
  } else if (CCTK_EQUALS(test_case, "Balsara4")) {
   
    CCTK_REAL rhol = 1.0;
    CCTK_REAL vxl = 0.999;
    CCTK_REAL vyl = 0.0;
    CCTK_REAL vzl = 0.0;
    CCTK_REAL pressl = 0.1;
    CCTK_REAL Bxl = 10.0;
    CCTK_REAL Byl = 7.0;
    CCTK_REAL Bzl = 7.0;

    CCTK_REAL rhor = 1.0;
    CCTK_REAL vxr = -0.999;
    CCTK_REAL vyr = 0.0;
    CCTK_REAL vzr = 0.0;
    CCTK_REAL pressr = 0.1;
    CCTK_REAL Bxr = 10.0;
    CCTK_REAL Byr = -7.0;
    CCTK_REAL Bzr = -7.0;    
	 
    if (CCTK_EQUALS(shock_dir, "x")) {
   
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = rhol;
            velx(p.I) = vxl;
            vely(p.I) = vyl;
            velz(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velx(p.I) = vxr;
            vely(p.I) = vyr;
            velz(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.x <= 0.0) { 
	       Avec_x(p.I) = Byl * (p.z) - Bzl * (p.y);  
	    } else { 
	       Avec_x(p.I) = Byr * (p.z) - Bzr * (p.y); } });


    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_z(p.I) = Bxr * (p.y);  });
    
    } else if (CCTK_EQUALS(shock_dir, "y")) {

    //x-->y, y-->z, z-->x
    Byl = 10.0;
    Bzl = 7.0;
    Bxl = 7.0;
    Byr = 10.0;
    Bzr = -7.0;
    Bxr = -7.0;

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.y <= 0.0) {
            rho(p.I) = rhol;
            vely(p.I) = vxl;
            velz(p.I) = vyl;
            velx(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            vely(p.I) = vxr;
            velz(p.I) = vyr;
            velx(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.y <= 0.0) {
               Avec_y(p.I) = Bzl * (p.x) - Bxl * (p.z);
            } else {
               Avec_y(p.I) = Bzr * (p.x) - Bxr * (p.z); } });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_x(p.I) = Byr * (p.z);  });

    } else if (CCTK_EQUALS(shock_dir, "z")) {
    
    //x-->z, y-->x, z-->y
    Bzl = 10.0;
    Bxl = 7.0;
    Byl = 7.0;
    Bzr = 10.0;
    Bxr = -7.0;
    Byr = -7.0;

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.z <= 0.0) {
            rho(p.I) = rhol;
            velz(p.I) = vxl;
            velx(p.I) = vyl;
            vely(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velz(p.I) = vxr;
            velx(p.I) = vyr;
            vely(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.z <= 0.0) {
               Avec_z(p.I) = Bxl * (p.y) - Byl * (p.x);
            } else {
               Avec_z(p.I) = Bxr * (p.y) - Byr * (p.x); } });


    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_y(p.I) = Bzr * (p.x);  });

    } else { CCTK_ERROR("Shock direction case not defined"); }
  
  } else if (CCTK_EQUALS(test_case, "Balsara5")) {
   
    CCTK_REAL rhol = 1.08;
    CCTK_REAL vxl = 0.4;
    CCTK_REAL vyl = 0.3;
    CCTK_REAL vzl = 0.2;
    CCTK_REAL pressl = 0.95;
    CCTK_REAL Bxl = 2.0;
    CCTK_REAL Byl = 0.3;
    CCTK_REAL Bzl = 0.3;

    CCTK_REAL rhor = 1.0;
    CCTK_REAL vxr = -0.45;
    CCTK_REAL vyr = -0.2;
    CCTK_REAL vzr = 0.2;
    CCTK_REAL pressr = 1.0;
    CCTK_REAL Bxr = 2.0;
    CCTK_REAL Byr = -0.7;
    CCTK_REAL Bzr = 0.5;    
	 
    if (CCTK_EQUALS(shock_dir, "x")) {
   
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = rhol;
            velx(p.I) = vxl;
            vely(p.I) = vyl;
            velz(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velx(p.I) = vxr;
            vely(p.I) = vyr;
            velz(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.x <= 0.0) { 
	       Avec_x(p.I) = Byl * (p.z) - Bzl * (p.y);  
	    } else { 
	       Avec_x(p.I) = Byr * (p.z) - Bzr * (p.y); } });


    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_z(p.I) = Bxr * (p.y);  });
    
    } else if (CCTK_EQUALS(shock_dir, "y")) {

    //x-->y, y-->z, z-->x
    Byl = 2.0;
    Bzl = 0.3;
    Bxl = 0.3;
    Byr = 2.0;
    Bzr = -0.7;
    Bxr = 0.5;

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.y <= 0.0) {
            rho(p.I) = rhol;
            vely(p.I) = vxl;
            velz(p.I) = vyl;
            velx(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            vely(p.I) = vxr;
            velz(p.I) = vyr;
            velx(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.y <= 0.0) {
               Avec_y(p.I) = Bzl * (p.x) - Bxl * (p.z);
            } else {
               Avec_y(p.I) = Bzr * (p.x) - Bxr * (p.z); } });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_x(p.I) = Byr * (p.z);  });

    } else if (CCTK_EQUALS(shock_dir, "z")) {
    
    //x-->z, y-->x, z-->y
    Bzl = 2.0;
    Bxl = 0.3;
    Byl = 0.3;
    Bzr = 2.0;
    Bxr = -0.7;
    Byr = 0.5;

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.z <= 0.0) {
            rho(p.I) = rhol;
            velz(p.I) = vxl;
            velx(p.I) = vyl;
            vely(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velz(p.I) = vxr;
            velx(p.I) = vyr;
            vely(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.z <= 0.0) {
               Avec_z(p.I) = Bxl * (p.y) - Byl * (p.x);
            } else {
               Avec_z(p.I) = Bxr * (p.y) - Byr * (p.x); } });


    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
               Avec_y(p.I) = Bzr * (p.x);  });

    } else { CCTK_ERROR("Shock direction case not defined"); }
    
  } else {
    CCTK_ERROR("Test case not defined");
  }
}

} // namespace AsterSeeds
