/* $Id: PDM_common_defines.h,v 1.2 2010/09/29 15:21:17 hspark Exp $ */
#ifndef PDM_COMMON_DEFINES_H
#define PDM_COMMON_DEFINES_H

/* Sequence of parameters is:
 * Ex
 * Ey
 * Ez
 * econv
 * alphatot
 * alpha1
 */

/* common variable defintion/mappings for auto-generated C code */

/* potential parameters */
double ex     	= params[ 0];
double ey     	= params[ 1];
double ez     	= params[ 2];
double econv  	= params[ 3];
double alphatot = params[ 4];
double alpha1   = params[ 5];

/* degrees of freedom */
double Xs1 = Xsi[0];
double Ys1 = Xsi[1];
double Zs1 = Xsi[2];

/* atom coordinates (reference) */
double X1 = Xa[0];	
double X2 = Xa[1];	
double X3 = Xa[2];	
double X4 = Xa[3];	
double X5 = Xa[4];	

double Y1 = Ya[0];	
double Y2 = Ya[1];	
double Y3 = Ya[2];	
double Y4 = Ya[3];	
double Y5 = Ya[4];	

double Z1 = Za[0];	
double Z2 = Za[1];	
double Z3 = Za[2];	
double Z4 = Za[3];	
double Z5 = Za[4];	

/* deformation */
double C11 = Cmat[0];
double C12 = Cmat[1];
double C13 = Cmat[2];
double C21 = Cmat[3];
double C22 = Cmat[4];
double C23 = Cmat[5];
double C31 = Cmat[6];
double C32 = Cmat[7];
double C33 = Cmat[8];

#endif /* PDM_COMMON_DEFINES_H */