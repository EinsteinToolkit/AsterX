#ifndef FM_IMPL_HXX
#define FM_IMPL_HXX

#include <cmath>
#include <cstdlib> // Needed for rand()
#include <cctk_Parameters.h>

namespace FMdisk
{

inline double GRHD_hm1(CCTK_REAL const& xcoord, CCTK_REAL const& ycoord, CCTK_REAL const& zcoord)
{

  DECLARE_CCTK_PARAMETERS;

  const double tmp_2 = 2*M*r_in;
  const double tmp_3 = ((a)*(a));
  const double tmp_4 = ((r_in)*(r_in)) + tmp_3;
  const double tmp_5 = -tmp_2 + tmp_4;
  const double tmp_6 = -tmp_3*tmp_5 + ((tmp_4)*(tmp_4));
  const double tmp_7 = ((r_at_max_density)*(r_at_max_density));
  const double tmp_8 = M*r_at_max_density;
  const double tmp_9 = a*sqrt(tmp_8);
  const double tmp_10 = tmp_7 - 3*tmp_8 + 2*tmp_9;
  const double tmp_11 = ((r_at_max_density)*(r_at_max_density)*(r_at_max_density)*(r_at_max_density)) + tmp_3*tmp_7 - 2*tmp_3*tmp_8 - tmp_9*(-tmp_3 + tmp_7);
  const double tmp_12 = M/((r_at_max_density)*(r_at_max_density)*(r_at_max_density));
  const double tmp_13 = 4*((tmp_11)*(tmp_11))*tmp_12/((tmp_10)*(tmp_10));
  const double tmp_14 = sqrt(((r_in)*(r_in)*(r_in)*(r_in))*tmp_13*tmp_5/((tmp_6)*(tmp_6)) + 1);
  const double tmp_16 = ((xcoord)*(xcoord)) + ((ycoord)*(ycoord)) + ((zcoord)*(zcoord));
  const double tmp_17 = 2*M*sqrt(tmp_16);
  const double tmp_19 = tmp_16 - tmp_17 + tmp_3;
  const double tmp_20 = ((zcoord)*(zcoord))/tmp_16;
  const double tmp_21 = tmp_16 + tmp_20*tmp_3;
  const double tmp_22 = 1 - tmp_20;
  const double tmp_23 = -tmp_19*tmp_22*tmp_3 + ((tmp_16 + tmp_3)*(tmp_16 + tmp_3));
  const double tmp_24 = sqrt(tmp_13*tmp_19*((tmp_21)*(tmp_21))/(tmp_22*((tmp_23)*(tmp_23))) + 1);
  const double tmp_25 = a*tmp_11*sqrt(tmp_12)/tmp_10;

  const double hm1 = -1 + sqrt(tmp_23*(tmp_24 + 1)/(tmp_19*tmp_21))*exp((1.0/2.0)*tmp_14 - tmp_17*tmp_25/tmp_23 + tmp_2*tmp_25/tmp_6 - 1.0/2.0*tmp_24)/sqrt(tmp_6*(tmp_14 + 1)/(((r_in)*(r_in))*tmp_5));

  return hm1;
}

void KerrSchild(CCTK_REAL const& xcoord, CCTK_REAL const& ycoord, CCTK_REAL const& zcoord,
                CCTK_REAL& alpha, CCTK_REAL& betaU0, CCTK_REAL& betaU1 , CCTK_REAL& betaU2,
                CCTK_REAL& gammaDD00, CCTK_REAL& gammaDD01, CCTK_REAL& gammaDD02, CCTK_REAL& gammaDD11, CCTK_REAL& gammaDD12, CCTK_REAL& gammaDD22,
                CCTK_REAL& KDD00, CCTK_REAL& KDD01, CCTK_REAL& KDD02, CCTK_REAL& KDD11, CCTK_REAL& KDD12, CCTK_REAL& KDD22)
{

  DECLARE_CCTK_PARAMETERS;

  const double FDPart3_0 = ((a)*(a));
  const double FDPart3_1 = ((zcoord)*(zcoord));
  const double FDPart3_2 = ((xcoord)*(xcoord));
  const double FDPart3_3 = ((ycoord)*(ycoord));
  const double FDPart3_4 = FDPart3_2 + FDPart3_3;
  const double FDPart3_5 = FDPart3_1 + FDPart3_4;
  const double FDPart3_6 = (1.0/(FDPart3_5));
  // FDPart3_7 is one on the z axis
  const double FDPart3_7 = FDPart3_1*FDPart3_6;
  const double FDPart3_8 = FDPart3_0*FDPart3_7 + FDPart3_5;
  const double FDPart3_9 = (1.0/(FDPart3_8));
  const double FDPart3_10 = sqrt(FDPart3_5);
  const double FDPart3_12 = 2*FDPart3_10*M;
  const double FDPart3_14 = FDPart3_12*FDPart3_9 + 1;
  const double FDPart3_15 = (1.0/(FDPart3_14));
  const double FDPart3_16 = M*xcoord;
  const double FDPart3_20 = 2*FDPart3_15*FDPart3_9;
  // FDPart3_22 is zero on the z axis
  const double FDPart3_22 = 1 - FDPart3_7;
  // FDPart3_24 is zero on the z axis
  const double FDPart3_24 = FDPart3_0*((FDPart3_14)*(FDPart3_14))*((FDPart3_22)*(FDPart3_22));
  const double FDPart3_25 = ((FDPart3_8)*(FDPart3_8));
  const double FDPart3_26 = FDPart3_0*FDPart3_12*FDPart3_22*FDPart3_9 + FDPart3_0 + FDPart3_5;
  // FDPart3_27 is zero on the z axis
  const double FDPart3_27 = FDPart3_14*FDPart3_22*FDPart3_26*FDPart3_8 - FDPart3_24*FDPart3_8;
  const double FDPart3_28 = (1.0/((FDPart3_27)*(FDPart3_27)*(FDPart3_27)));
  const double FDPart3_29 = FDPart3_14*FDPart3_22*FDPart3_26 - FDPart3_24;
  const double FDPart3_30 = (1.0/(FDPart3_14*FDPart3_22*FDPart3_25*FDPart3_26*FDPart3_28*FDPart3_29 - FDPart3_24*FDPart3_25*FDPart3_28*FDPart3_29));
  const double FDPart3_31 = (1.0/((FDPart3_27)*(FDPart3_27)));
  const double FDPart3_32 = FDPart3_29*FDPart3_30*FDPart3_31*FDPart3_8;
  const double FDPart3_33 = FDPart3_14*FDPart3_32;
  const double FDPart3_34 = xcoord*ycoord;
  const double FDPart3_35 = (1.0/(FDPart3_4));
  const double FDPart3_36 = FDPart3_22*FDPart3_35;
  const double FDPart3_37 = (1.0/(FDPart3_10));
  const double FDPart3_39 = -FDPart3_14*FDPart3_37*a;
  const double FDPart3_41 = 2*FDPart3_32*FDPart3_34*FDPart3_36*FDPart3_39;
  const double FDPart3_43 = (1.0/((FDPart3_4)*(FDPart3_4)));
  const double FDPart3_44 = FDPart3_22*FDPart3_26*FDPart3_43;
  const double FDPart3_46 = (1.0/(FDPart3_22));
  const double FDPart3_47 = FDPart3_30*(FDPart3_14*FDPart3_22*FDPart3_25*FDPart3_26*FDPart3_31 - FDPart3_24*FDPart3_25*FDPart3_31);
  const double FDPart3_49 = FDPart3_46*FDPart3_47/((FDPart3_5)*(FDPart3_5)*(FDPart3_5));
  const double FDPart3_50 = FDPart3_33*FDPart3_6;
  const double FDPart3_54 = FDPart3_37*zcoord;
  const double FDPart3_56 = -FDPart3_14*FDPart3_32*FDPart3_36*FDPart3_54*a;
  const double FDPart3_57 = pow(FDPart3_5, -3.0/2.0);
  const double FDPart3_58 = FDPart3_1*FDPart3_57;
  const double FDPart3_59 = FDPart3_37 - FDPart3_58;
  const double FDPart3_60 = FDPart3_46*FDPart3_47*FDPart3_57*FDPart3_59;
  const double FDPart3_62 = FDPart3_46*((FDPart3_59)*(FDPart3_59));
  const double FDPart3_63 = (1.0/((FDPart3_5)*(FDPart3_5)));
  const double FDPart3_64 = FDPart3_1*FDPart3_2*FDPart3_63;
  const double FDPart3_65 = sqrt(FDPart3_14);
  const double FDPart3_67 = acos(FDPart3_54);
  const double FDPart3_68 = cos(2*FDPart3_67);
  const double FDPart3_70 = 2*FDPart3_2;
  const double FDPart3_71 = 2*FDPart3_3;
  const double FDPart3_72 = 2*FDPart3_1 + FDPart3_70 + FDPart3_71;
  const double FDPart3_73 = FDPart3_0*FDPart3_68 + FDPart3_0 + FDPart3_72;
  const double FDPart3_74 = (1.0/(4*FDPart3_10*M + FDPart3_73));
  const double FDPart3_75 = FDPart3_65*FDPart3_74;
  const double FDPart3_76 = FDPart3_75*M;
  const double FDPart3_78 = 4*FDPart3_46*FDPart3_76;
  const double FDPart3_79 = (1.0/(FDPart3_73));
  const double FDPart3_81 = 16*FDPart3_0*FDPart3_79;
  const double FDPart3_82 = FDPart3_76*FDPart3_81;
  const double FDPart3_83 = (1.0/((FDPart3_73)*(FDPart3_73)));
  const double FDPart3_84 = FDPart3_0*FDPart3_68 + FDPart3_0 - FDPart3_72;
  const double FDPart3_86 = FDPart3_36*FDPart3_65*FDPart3_83*FDPart3_84;
  const double FDPart3_87 = FDPart3_37*FDPart3_86*a;
  const double FDPart3_88 = 4*FDPart3_16*FDPart3_87*ycoord;
  const double FDPart3_89 = ((a)*(a)*(a));
  const double FDPart3_90 = FDPart3_16*FDPart3_75;
  const double FDPart3_91 = FDPart3_90*ycoord;
  const double FDPart3_92 = 16*FDPart3_36*FDPart3_58*FDPart3_79*FDPart3_89*FDPart3_91;
  const double FDPart3_94 = FDPart3_83*FDPart3_84*(FDPart3_12 + FDPart3_73);
  const double FDPart3_95 = 4*FDPart3_76*FDPart3_94;
  const double FDPart3_96 = ((a)*(a)*(a)*(a));
  const double FDPart3_98 = FDPart3_22*FDPart3_43*FDPart3_83*(4*FDPart3_0*FDPart3_10*FDPart3_68*(FDPart3_0 + FDPart3_10*(2*FDPart3_10 + M)) + 4*FDPart3_0*FDPart3_5*(2*FDPart3_10 - M) + 8*pow(FDPart3_5, 5.0/2.0) + FDPart3_96*(FDPart3_10 - M)*cos(4*FDPart3_67) + FDPart3_96*(3*FDPart3_10 + M));
  const double FDPart3_99 = FDPart3_10*FDPart3_76*FDPart3_98;
  const double FDPart3_102 = FDPart3_1*FDPart3_63*FDPart3_91;
  const double FDPart3_104 = 8*FDPart3_79*FDPart3_89;
  const double FDPart3_105 = FDPart3_104*FDPart3_58*FDPart3_76;
  const double FDPart3_106 = 4*FDPart3_6*FDPart3_94;
  const double FDPart3_107 = ((zcoord)*(zcoord)*(zcoord));
  const double FDPart3_108 = FDPart3_54*FDPart3_59;
  const double FDPart3_109 = 4*FDPart3_108*FDPart3_46;
  const double FDPart3_110 = 8*FDPart3_0*FDPart3_79;
  const double FDPart3_111 = FDPart3_75*M*ycoord;
  const double FDPart3_112 = FDPart3_104*FDPart3_36*FDPart3_59*zcoord;
  const double FDPart3_113 = FDPart3_1*FDPart3_3*FDPart3_63;

  alpha = sqrt(FDPart3_15);
  betaU0 = 2*FDPart3_15*FDPart3_16*FDPart3_9;
  betaU1 = FDPart3_20*M*ycoord;
  betaU2 = FDPart3_20*M*zcoord;

  gammaDD00 = FDPart3_1*FDPart3_2*FDPart3_49 + FDPart3_2*FDPart3_33*FDPart3_6 + FDPart3_3*FDPart3_32*FDPart3_44 - FDPart3_41;
  gammaDD01 = FDPart3_1*FDPart3_34*FDPart3_49 + FDPart3_2*FDPart3_32*FDPart3_36*FDPart3_39 - FDPart3_3*FDPart3_32*FDPart3_36*FDPart3_39 - FDPart3_32*FDPart3_34*FDPart3_44 + FDPart3_34*FDPart3_50;
  gammaDD02 = FDPart3_14*FDPart3_29*FDPart3_30*FDPart3_31*FDPart3_6*FDPart3_8*xcoord*zcoord - FDPart3_56*ycoord - FDPart3_60*xcoord*zcoord;
  gammaDD11 = FDPart3_1*FDPart3_3*FDPart3_49 + FDPart3_2*FDPart3_32*FDPart3_44 + FDPart3_3*FDPart3_50 + FDPart3_41;
  gammaDD12 = FDPart3_50*ycoord*zcoord + FDPart3_56*xcoord - FDPart3_60*ycoord*zcoord;
  gammaDD22 = FDPart3_33*FDPart3_7 + FDPart3_47*FDPart3_62;

  KDD00 = FDPart3_2*FDPart3_6*FDPart3_95 + FDPart3_64*FDPart3_78 + FDPart3_64*FDPart3_82 + FDPart3_71*FDPart3_99 + FDPart3_88 + FDPart3_92;
  KDD01 = 4*FDPart3_102*FDPart3_46 + FDPart3_102*FDPart3_81 - FDPart3_105*FDPart3_2*FDPart3_36 + FDPart3_105*FDPart3_3*FDPart3_36 + FDPart3_106*FDPart3_91 - FDPart3_12*FDPart3_34*FDPart3_75*FDPart3_98 - FDPart3_70*FDPart3_87*M + FDPart3_71*FDPart3_87*M;
  KDD02 = 8*FDPart3_0*FDPart3_107*FDPart3_63*FDPart3_65*FDPart3_74*FDPart3_79*M*xcoord - FDPart3_108*FDPart3_110*FDPart3_90 - FDPart3_109*FDPart3_90 - FDPart3_111*FDPart3_112 + 2*FDPart3_22*FDPart3_35*FDPart3_37*FDPart3_65*FDPart3_83*FDPart3_84*M*a*ycoord*zcoord + 4*FDPart3_6*FDPart3_65*FDPart3_74*FDPart3_83*FDPart3_84*M*xcoord*zcoord*(FDPart3_12 + FDPart3_73);
  KDD11 = FDPart3_113*FDPart3_78 + FDPart3_113*FDPart3_82 + FDPart3_3*FDPart3_6*FDPart3_95 + FDPart3_70*FDPart3_99 - FDPart3_88 - FDPart3_92;
  KDD12 = FDPart3_106*FDPart3_111*zcoord + FDPart3_107*FDPart3_110*FDPart3_111*FDPart3_63 - FDPart3_108*FDPart3_110*FDPart3_111 - FDPart3_109*FDPart3_111 + FDPart3_112*FDPart3_90 - 2*FDPart3_16*FDPart3_54*FDPart3_86*a;
  KDD22 = -FDPart3_1*FDPart3_37*FDPart3_59*FDPart3_82 + 4*FDPart3_5*FDPart3_62*FDPart3_76 + FDPart3_7*FDPart3_95;

}

void GRHD_velocities(CCTK_REAL const& xcoord, CCTK_REAL const& ycoord, CCTK_REAL const& zcoord,
                     CCTK_REAL &Valencia3velocityU0GF, CCTK_REAL &Valencia3velocityU1GF, CCTK_REAL &Valencia3velocityU2GF)
{

  DECLARE_CCTK_PARAMETERS;

  const double FDPart3_0 = ((a)*(a));
  const double FDPart3_2 = ((xcoord)*(xcoord)) + ((ycoord)*(ycoord)) + ((zcoord)*(zcoord));
  const double FDPart3_3 = ((zcoord)*(zcoord))/FDPart3_2;
  const double FDPart3_4 = FDPart3_0*FDPart3_3 + FDPart3_2;
  const double FDPart3_5 = (1.0/(FDPart3_4));
  const double FDPart3_7 = 2*sqrt(FDPart3_2)*M;
  const double FDPart3_9 = (1.0/(FDPart3_5*FDPart3_7 + 1));
  const double FDPart3_10 = (1.0/sqrt(FDPart3_9));
  const double FDPart3_12 = 2*FDPart3_5*FDPart3_9*M;
  const double FDPart3_15 = 1 - FDPart3_3;
  const double FDPart3_16 = FDPart3_0 + FDPart3_2 - FDPart3_7;
  const double FDPart3_17 = FDPart3_0*FDPart3_15*FDPart3_16 - ((FDPart3_0 + FDPart3_2)*(FDPart3_0 + FDPart3_2));
  const double FDPart3_19 = sqrt(-FDPart3_17*FDPart3_5);
  const double FDPart3_20 = ((r_at_max_density)*(r_at_max_density));
  const double FDPart3_21 = M*r_at_max_density;
  const double FDPart3_22 = sqrt(FDPart3_21)*a;
  const double FDPart3_23 = (1.0/(FDPart3_15));
  const double FDPart3_24 = (1.0/2.0)*sqrt(4*FDPart3_16*FDPart3_23*((FDPart3_4)*(FDPart3_4))*M*((FDPart3_0*FDPart3_20 - 2*FDPart3_0*FDPart3_21 - FDPart3_22*(-FDPart3_0 + FDPart3_20) + ((r_at_max_density)*(r_at_max_density)*(r_at_max_density)*(r_at_max_density)))*(FDPart3_0*FDPart3_20 - 2*FDPart3_0*FDPart3_21 - FDPart3_22*(-FDPart3_0 + FDPart3_20) + ((r_at_max_density)*(r_at_max_density)*(r_at_max_density)*(r_at_max_density))))/(((FDPart3_17)*(FDPart3_17))*((r_at_max_density)*(r_at_max_density)*(r_at_max_density))*((FDPart3_20 - 3*FDPart3_21 + 2*FDPart3_22)*(FDPart3_20 - 3*FDPart3_21 + 2*FDPart3_22))) + 1);
  const double FDPart3_25 = sqrt(FDPart3_24 - 1.0/2.0);
  const double FDPart3_26 = fabs(sqrt(FDPart3_15));
  const double FDPart3_27 = FDPart3_19*FDPart3_25*FDPart3_26;
  const double FDPart3_28 = (1.0/(FDPart3_16));
  const double FDPart3_29 = FDPart3_28*FDPart3_5*FDPart3_7*a;
  const double FDPart3_30 = -1/FDPart3_17;
  const double FDPart3_31 = -FDPart3_27*FDPart3_30*FDPart3_7*a - sqrt(FDPart3_4)*sqrt(FDPart3_16*FDPart3_30)*sqrt(FDPart3_24 + 1.0/2.0);
  const double FDPart3_32 = (FDPart3_19*FDPart3_25*FDPart3_26*(-4*FDPart3_0*FDPart3_2*FDPart3_28*FDPart3_30*FDPart3_5*((M)*(M)) + FDPart3_23*FDPart3_30*FDPart3_4) - FDPart3_29*FDPart3_31)/(FDPart3_17*FDPart3_28*FDPart3_31*FDPart3_5 - FDPart3_27*FDPart3_29);

  Valencia3velocityU0GF = FDPart3_10*FDPart3_12*xcoord - FDPart3_10*FDPart3_32*ycoord;
  Valencia3velocityU1GF = FDPart3_10*FDPart3_12*ycoord + FDPart3_10*FDPart3_32*xcoord;
  Valencia3velocityU2GF = FDPart3_10*FDPart3_12*zcoord;
}

void GRHD_perturb_pressure(double &press, double &eps, double const &rho) {

  DECLARE_CCTK_PARAMETERS;

  // Generate random number in range [0,1),
  // snippet courtesy http://daviddeley.com/random/crandom.htm
  const double random_number_between_0_and_1 = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );

  const double random_number_between_min_and_max = random_min + (random_max - random_min)*random_number_between_0_and_1;

  press = press*(1.0 + random_number_between_min_and_max);

  // Add 1e-300 to rho to avoid division by zero when density is zero.
  eps = press / ((rho + 1e-300) * (gamma - 1.0));

}

//void GRMHD_set_A(double const &press, double const &xtilde, double const &ytilde, double &Ax, double &Ay, double &Az) {
//
//  DECLARE_CCTK_PARAMETERS;
//
//  Ax = - A_b * pow(fmax(press-press_cut,0.),A_n) * ytilde;
//  Ay =   A_b * pow(fmax(press-press_cut,0.),A_n) * xtilde;
//  Az = 0.;
//}

void GRMHD_set_A(double const &press, double const &rho, double const &xtilde, double const &ytilde, double &Ax, double &Ay, double &Az) {

  DECLARE_CCTK_PARAMETERS;

  const double rcyl = fmax(sqrt(xtilde*xtilde+ytilde*ytilde),1e-15);

  // If A_{\theta} = A_r = 0:
  // A_x = (-y/rcyl/rcyl)*A_{\phi}
  // A_y = (x/rcyl/rcyl)*A_{\phi}
  // A_z = 0
  
  if (use_pressure) {

   Ax = - A_b * pow(fmax(press-press_cut,0.),A_n) * ytilde * pow(rcyl,A_c-2.);
   Ay =   A_b * pow(fmax(press-press_cut,0.),A_n) * xtilde * pow(rcyl,A_c-2.);
   Az = 0.;

  } else {

   Ax = - A_b * pow(fmax(rho-rho_cut,0.),A_n) * ytilde * pow(rcyl,A_c-2.);
   Ay =   A_b * pow(fmax(rho-rho_cut,0.),A_n) * xtilde * pow(rcyl,A_c-2.);
   Az = 0.;

  }

}

} // end namespace FMdisk

#endif // FM_IMPL_HXX
