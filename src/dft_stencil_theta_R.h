/* This file was automatically generated.  Do not edit! */
double StenTheta_R_GetWeightFromSten(double rsq,double R);
int StenTheta_R_NquadPtsGauss(double r);
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAS_VALUES_H)
#include <values.h>
#include <unistd.h>
#include <string.h>
#endif
#include "mpi.h"
#include "az_aztec.h"
#include "az_aztec_defs.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
extern int Ndim;
int StenTheta_R_NquadPtsBoundary();
int StenTheta_R_Njcomp();
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
#define PI    3.141592653589793238462643383279502884197169399375
double StenTheta_R_sten_vol(int i);
#define NCOMP_MAX 5
extern double HS_diam[NCOMP_MAX];
double StenTheta_R_sten_rad(int icomp);
