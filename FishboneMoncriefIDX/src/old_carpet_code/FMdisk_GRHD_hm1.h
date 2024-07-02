/*
 *  Original SymPy expression:
 *  "hm1 = -1 + sqrt((-a**2*(-zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + 1)*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2) + (a**2 + xcoord**2 + ycoord**2 + zcoord**2)**2)*(sqrt(4*M*(-2*M*a**2*r_at_max_density + a**2*r_at_max_density**2 - a*sqrt(M*r_at_max_density)*(-a**2 + r_at_max_density**2) + r_at_max_density**4)**2*(a**2*zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + xcoord**2 + ycoord**2 + zcoord**2)**2*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2)/(r_at_max_density**3*(-zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + 1)*(-a**2*(-zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + 1)*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2) + (a**2 + xcoord**2 + ycoord**2 + zcoord**2)**2)**2*(-3*M*r_at_max_density + 2*a*sqrt(M*r_at_max_density) + r_at_max_density**2)**2) + 1) + 1)/((a**2*zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + xcoord**2 + ycoord**2 + zcoord**2)*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2)))*exp(2*M*a*r_in*sqrt(M/r_at_max_density**3)*(-2*M*a**2*r_at_max_density + a**2*r_at_max_density**2 - a*sqrt(M*r_at_max_density)*(-a**2 + r_at_max_density**2) + r_at_max_density**4)/((-a**2*(-2*M*r_in + a**2 + r_in**2) + (a**2 + r_in**2)**2)*(-3*M*r_at_max_density + 2*a*sqrt(M*r_at_max_density) + r_at_max_density**2)) - 2*M*a*sqrt(M/r_at_max_density**3)*sqrt(xcoord**2 + ycoord**2 + zcoord**2)*(-2*M*a**2*r_at_max_density + a**2*r_at_max_density**2 - a*sqrt(M*r_at_max_density)*(-a**2 + r_at_max_density**2) + r_at_max_density**4)/((-a**2*(-zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + 1)*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2) + (a**2 + xcoord**2 + ycoord**2 + zcoord**2)**2)*(-3*M*r_at_max_density + 2*a*sqrt(M*r_at_max_density) + r_at_max_density**2)) + sqrt(4*M*r_in**4*(-2*M*r_in + a**2 + r_in**2)*(-2*M*a**2*r_at_max_density + a**2*r_at_max_density**2 - a*sqrt(M*r_at_max_density)*(-a**2 + r_at_max_density**2) + r_at_max_density**4)**2/(r_at_max_density**3*(-a**2*(-2*M*r_in + a**2 + r_in**2) + (a**2 + r_in**2)**2)**2*(-3*M*r_at_max_density + 2*a*sqrt(M*r_at_max_density) + r_at_max_density**2)**2) + 1)/2 - sqrt(4*M*(-2*M*a**2*r_at_max_density + a**2*r_at_max_density**2 - a*sqrt(M*r_at_max_density)*(-a**2 + r_at_max_density**2) + r_at_max_density**4)**2*(a**2*zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + xcoord**2 + ycoord**2 + zcoord**2)**2*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2)/(r_at_max_density**3*(-zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + 1)*(-a**2*(-zcoord**2/(xcoord**2 + ycoord**2 + zcoord**2) + 1)*(-2*M*sqrt(xcoord**2 + ycoord**2 + zcoord**2) + a**2 + xcoord**2 + ycoord**2 + zcoord**2) + (a**2 + xcoord**2 + ycoord**2 + zcoord**2)**2)**2*(-3*M*r_at_max_density + 2*a*sqrt(M*r_at_max_density) + r_at_max_density**2)**2) + 1)/2)/sqrt((-a**2*(-2*M*r_in + a**2 + r_in**2) + (a**2 + r_in**2)**2)*(sqrt(4*M*r_in**4*(-2*M*r_in + a**2 + r_in**2)*(-2*M*a**2*r_at_max_density + a**2*r_at_max_density**2 - a*sqrt(M*r_at_max_density)*(-a**2 + r_at_max_density**2) + r_at_max_density**4)**2/(r_at_max_density**3*(-a**2*(-2*M*r_in + a**2 + r_in**2) + (a**2 + r_in**2)**2)**2*(-3*M*r_at_max_density + 2*a*sqrt(M*r_at_max_density) + r_at_max_density**2)**2) + 1) + 1)/(r_in**2*(-2*M*r_in + a**2 + r_in**2)))"
 */
{
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
  hm1 = -1 + sqrt(tmp_23*(tmp_24 + 1)/(tmp_19*tmp_21))*exp((1.0/2.0)*tmp_14 - tmp_17*tmp_25/tmp_23 + tmp_2*tmp_25/tmp_6 - 1.0/2.0*tmp_24)/sqrt(tmp_6*(tmp_14 + 1)/(((r_in)*(r_in))*tmp_5));
}
