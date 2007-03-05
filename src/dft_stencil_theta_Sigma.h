/* This file was automatically generated.  Do not edit! */
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
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define PI    M_PI
double StenTheta_Sigma_GetWeightFromSten(double rsq,double R);
int StenTheta_Sigma_NquadPtsGauss(double r);
extern int Ndim;
int StenTheta_Sigma_NquadPtsBoundary();
int StenTheta_Sigma_Njcomp();
double StenTheta_Sigma_sten_vol(int icomp);
#define NCOMP_MAX 5
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
double StenTheta_Sigma_sten_rad(int icomp);
