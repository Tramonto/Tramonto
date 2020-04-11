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
#include "az_aztec_defs.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define PI    3.141592653589793238462643383279502884197169399375
double deltaC_MSA(double r,int i,int j);
double StenTheta_RPMmsa_GetWeightFromSten(int icomp,int jcomp,double rsq,int ngpu,double *gpu,double *gwu);
int StenTheta_RPMmsa_NquadPtsGaussIntegrand(double r);
int StenTheta_RPMmsa_NquadPtsGauss(double r);
extern int Ndim;
int StenTheta_RPMmsa_NquadPtsBoundary();
extern int Ncomp;
int StenTheta_RPMmsa_Njcomp();
double deltaC_MSA_int(double r,int i,int j);
double StenTheta_RPMmsa_sten_vol(int i,int j);
#define NCOMP_MAX 6
extern double HS_diam[NCOMP_MAX];
double StenTheta_RPMmsa_sten_rad(int icomp);
