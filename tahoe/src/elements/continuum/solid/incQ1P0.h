#ifndef INCQ1P0_H
#define INCQ1P0_H

#ifdef __cplusplus
extern "C" {
#endif

/* Sequence of parameters is:
 * epsilon
*/

/* function to compute first derivative of free energy wrt to the stretch tensor C */
void mech_pk2_q1p0(const double epsilon, const double* Xsi, const double* Cmat, const double* Fmat, double J, double I1, double* dUdCmech);
void me_pk2_q1p0(const double epsilon, const double* Xsi, const double* Cmat, const double* Fmat, double J, double* dUdCmechelec);

#ifdef __cplusplus
}
#endif

#endif /* INCQ1P0_H */
