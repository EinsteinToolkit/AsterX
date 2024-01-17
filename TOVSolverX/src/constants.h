/* From EinsteinBase/Constants/src/constants.h */
/* $Header$ */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifdef CCODE
#define CONSTANT_c_SI  2.99792458e8
#define CONSTANT_c_cgi 2.99792458e10
#define CONSTANT_c_C   1.0e0
#endif
#ifdef FCODE
#define CONSTANT_c_SI  2.99792458d8
#define CONSTANT_c_cgi 2.9979245800d10
#define CONSTANT_c_C   1.0d0
#endif

#ifdef CCODE
#define CONSTANT_G_SI  6.6732e-11
#define CONSTANT_G_cgi 6.6732e-8
#define CONSTANT_G_C   1.0e0
#endif
#ifdef FCODE
#define CONSTANT_G_SI  6.6732d-11
#define CONSTANT_G_cgi 6.6732d-8
#define CONSTANT_G_C   1.0e0
#endif

#ifdef CCODE
#define CONSTANT_Msolar_SI  1.987e30
#define CONSTANT_Msolar_cgi 1.987e33
#define CONSTANT_Msolar_C   1.0d0
#endif
#ifdef FCODE
#define CONSTANT_Msolar_SI  1.987d30
#define CONSTANT_Msolar_cgi 1.987d33
#define CONSTANT_Msolar_C   1.0d0
#endif

#define CONSTANT

#endif
