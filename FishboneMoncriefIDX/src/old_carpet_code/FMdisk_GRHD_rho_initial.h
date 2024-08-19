{
  /*
   * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
   */
  /*
   *  Original SymPy expression:
   */
  const double xcoord = xcoordGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
  const double ycoord = ycoordGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
  const double zcoord = zcoordGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
  /*
   * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
   */
  /*
   *  Original SymPy expression:
   *  "rho_initialGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = exp(log((-1 + sqrt((-a**2*(-zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + 1)*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2) + (a**2 + xcoord**2 + ycoord**2 + zcoord**2)**2)*(sqrt(4*M*(-2*M*a**2*r_at_max_density + a**2*r_at_max_density**2 - a*sqrt(M*r_at_max_density)*(-a**2 + r_at_max_density**2) + r_at_max_density**4)**2*(a**2*zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + xcoord**2 + ycoord**2 + zcoord**2)**2*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2)/(r_at_max_density**3*(-zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + 1)*(-a**2*(-zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + 1)*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2) + (a**2 + xcoord**2 + ycoord**2 + zcoord**2)**2)**2*(-3*M*r_at_max_density + 2*a*sqrt(M*r_at_max_density) + r_at_max_density**2)**2) + 1) + 1)/((a**2*zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + xcoord**2 + ycoord**2 + zcoord**2)*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2)))*exp(2*M*a*r_in*sqrt(M/r_at_max_density**3)*(-2*M*a**2*r_at_max_density + a**2*r_at_max_density**2 - a*sqrt(M*r_at_max_density)*(-a**2 + r_at_max_density**2) + r_at_max_density**4)/((-a**2*(-2*M*r_in + a**2 + r_in**2) + (a**2 + r_in**2)**2)*(-3*M*r_at_max_density + 2*a*sqrt(M*r_at_max_density) + r_at_max_density**2)) - 2*M*a*sqrt(M/r_at_max_density**3)*sqrt(xcoord**2 + ycoord**2 + zcoord**2)*(-2*M*a**2*r_at_max_density + a**2*r_at_max_density**2 - a*sqrt(M*r_at_max_density)*(-a**2 + r_at_max_density**2) + r_at_max_density**4)/((-a**2*(-zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + 1)*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2) + (a**2 + xcoord**2 + ycoord**2 + zcoord**2)**2)*(-3*M*r_at_max_density + 2*a*sqrt(M*r_at_max_density) + r_at_max_density**2)) + sqrt(4*M*r_in**4*(-2*M*r_in + a**2 + r_in**2)*(-2*M*a**2*r_at_max_density + a**2*r_at_max_density**2 - a*sqrt(M*r_at_max_density)*(-a**2 + r_at_max_density**2) + r_at_max_density**4)**2/(r_at_max_density**3*(-a**2*(-2*M*r_in + a**2 + r_in**2) + (a**2 + r_in**2)**2)**2*(-3*M*r_at_max_density + 2*a*sqrt(M*r_at_max_density) + r_at_max_density**2)**2) + 1)/2 - sqrt(4*M*(-2*M*a**2*r_at_max_density + a**2*r_at_max_density**2 - a*sqrt(M*r_at_max_density)*(-a**2 + r_at_max_density**2) + r_at_max_density**4)**2*(a**2*zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + xcoord**2 + ycoord**2 + zcoord**2)**2*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2)/(r_at_max_density**3*(-zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + 1)*(-a**2*(-zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + 1)*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2) + (a**2 + xcoord**2 + ycoord**2 + zcoord**2)**2)**2*(-3*M*r_at_max_density + 2*a*sqrt(M*r_at_max_density) + r_at_max_density**2)**2) + 1)/2)/sqrt((-a**2*(-2*M*r_in + a**2 + r_in**2) + (a**2 + r_in**2)**2)*(sqrt(4*M*r_in**4*(-2*M*r_in + a**2 + r_in**2)*(-2*M*a**2*r_at_max_density + a**2*r_at_max_density**2 - a*sqrt(M*r_at_max_density)*(-a**2 + r_at_max_density**2) + r_at_max_density**4)**2/(r_at_max_density**3*(-a**2*(-2*M*r_in + a**2 + r_in**2) + (a**2 + r_in**2)**2)**2*(-3*M*r_at_max_density + 2*a*sqrt(M*r_at_max_density) + r_at_max_density**2)**2) + 1) + 1)/(r_in**2*(-2*M*r_in + a**2 + r_in**2))))*(gamma - 1)/(gamma*kappa))/(gamma - 1))"
   */
  const double FDPart3_3 = 2*M*r_in;
  const double FDPart3_4 = ((a)*(a));
  const double FDPart3_5 = FDPart3_4 + ((r_in)*(r_in));
  const double FDPart3_6 = -FDPart3_3 + FDPart3_5;
  const double FDPart3_7 = -FDPart3_4*FDPart3_6 + ((FDPart3_5)*(FDPart3_5));
  const double FDPart3_8 = ((r_at_max_density)*(r_at_max_density));
  const double FDPart3_9 = M*r_at_max_density;
  const double FDPart3_10 = sqrt(FDPart3_9)*a;
  const double FDPart3_11 = 2*FDPart3_10 + FDPart3_8 - 3*FDPart3_9;
  const double FDPart3_12 = -FDPart3_10*(-FDPart3_4 + FDPart3_8) + FDPart3_4*FDPart3_8 - 2*FDPart3_4*FDPart3_9 + ((r_at_max_density)*(r_at_max_density)*(r_at_max_density)*(r_at_max_density));
  const double FDPart3_13 = M/((r_at_max_density)*(r_at_max_density)*(r_at_max_density));
  const double FDPart3_14 = 4*((FDPart3_12)*(FDPart3_12))*FDPart3_13/((FDPart3_11)*(FDPart3_11));
  const double FDPart3_15 = sqrt(FDPart3_14*FDPart3_6*((r_in)*(r_in)*(r_in)*(r_in))/((FDPart3_7)*(FDPart3_7)) + 1);
  const double FDPart3_17 = ((xcoord)*(xcoord)) + ((ycoord)*(ycoord)) + ((zcoord)*(zcoord));
  const double FDPart3_18 = 2*sqrt(FDPart3_17)*M;
  const double FDPart3_20 = FDPart3_17 - FDPart3_18 + FDPart3_4;
  const double FDPart3_21 = ((zcoord)*(zcoord))/FDPart3_17;
  const double FDPart3_22 = FDPart3_17 + FDPart3_21*FDPart3_4;
  const double FDPart3_23 = 1 - FDPart3_21;
  const double FDPart3_24 = -FDPart3_20*FDPart3_23*FDPart3_4 + ((FDPart3_17 + FDPart3_4)*(FDPart3_17 + FDPart3_4));
  const double FDPart3_25 = sqrt(FDPart3_14*FDPart3_20*((FDPart3_22)*(FDPart3_22))/(FDPart3_23*((FDPart3_24)*(FDPart3_24))) + 1);
  const double FDPart3_26 = FDPart3_12*sqrt(FDPart3_13)*a/FDPart3_11;
  rho_initialGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = exp(log((gamma - 1)*(sqrt(FDPart3_24*(FDPart3_25 + 1)/(FDPart3_20*FDPart3_22))*exp((1.0/2.0)*FDPart3_15 - FDPart3_18*FDPart3_26/FDPart3_24 - 1.0/2.0*FDPart3_25 + FDPart3_26*FDPart3_3/FDPart3_7)/sqrt(FDPart3_7*(FDPart3_15 + 1)/(FDPart3_6*((r_in)*(r_in)))) - 1)/(gamma*kappa))/(gamma - 1));
}