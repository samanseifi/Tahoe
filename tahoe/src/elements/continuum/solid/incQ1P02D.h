#ifndef INCQ1P02D_H
#define INCQ1P02D_H

#ifdef __cplusplus
extern "C" {
#endif

/* Sequence of parameters is:
 * epsilon
 * mu
 * Nrig
 * lambda
 */

/* function to compute first derivative of free energy wrt to the stretch tensor C */
void mech_pk2_q1p02D(const double* params, const double* Xsi, const double* Cmat, const double* Fmat, double J, double I1, double* dUdCmech);
void me_pk2_q1p02D(const double* params, const double* Xsi, const double* Cmat, const double* Fmat, double J, double* dUdCmechelec);

#ifdef __cplusplus
}
#endif

#endif /* INCQ1P02D_H */
