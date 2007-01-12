/* This file was automatically generated.  Do not edit! */
double StenDelta_R_GetWeightFromSten(double rsq,double R);
int StenDelta_R_NquadPtsGauss(double r);
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAS_VALUES_H)
#include <values.h>
#include <unistd.h>
#endif
#include "mpi.h"
#include "az_aztec.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
extern int Ndim;
int StenDelta_R_NquadPtsBoundary();
int StenDelta_R_Njcomp();
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
#define PI    M_PI
double StenDelta_R_sten_vol(int icomp);
#define NCOMP_MAX 5
extern double HS_diam[NCOMP_MAX];
double StenDelta_R_sten_rad(int icomp);
