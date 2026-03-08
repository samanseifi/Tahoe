#ifndef INCQ1P0_ELEC_H
#define INCQ1P0_ELEC_H

#ifdef __cplusplus
extern "C" {
#endif

/* Sequence of parameters is:
 * epsilon
*/
 /* function to get electric displacement */
void elec_pk2_q1p0(const double epsilon, const double *Xsi, const double* Cmat, const double* Fmat, double J, double* dUdCelec); 

#ifdef __cplusplus
}
#endif

#endif /* INCQ1P0_ELEC_H */
