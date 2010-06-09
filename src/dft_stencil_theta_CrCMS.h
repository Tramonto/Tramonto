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
#define PI    3.141592653589793238462643383279502884197169399375
double int_cr(double r_low,double r_upp,double slope_dr,int icomp,int jcomp,int irmin,double zsq,double *rx_low);
extern double ***Rism_cr;
extern int Last_nz_cr;
extern double Deltar_cr;
double StenTheta_CrCMS_GetWeightFromSten(int icomp,int jcomp,double rsq,double R);
int StenTheta_CrCMS_NquadPtsGauss(double r);
extern int Ndim;
int StenTheta_CrCMS_NquadPtsBoundary();
extern int Ncomp;
int StenTheta_CrCMS_Njcomp();
#define NO_RENORMALIZATION_FLAG -888
double StenTheta_CrCMS_sten_vol(int icomp,int jcomp);
#define NCOMP_MAX 5
extern double Cr_rad[NCOMP_MAX][NCOMP_MAX];
double StenTheta_CrCMS_sten_rad(int icomp,int jcomp);
