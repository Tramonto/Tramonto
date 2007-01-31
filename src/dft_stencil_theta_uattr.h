/* This file was automatically generated.  Do not edit! */
double pairPot_ATT_CS_switch(double r,int icomp,int jcomp,int typePairPot);
double StenTheta_uattr_GetWeightFromSten(int icomp,int jcomp,double rsq,int ngpu,double *gpu,double *gwu);
int StenTheta_uattr_NquadPtsGaussIntegrand(double r);
int StenTheta_uattr_NquadPtsGauss(double r);
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
int StenTheta_uattr_NquadPtsBoundary();
extern int Ncomp;
int StenTheta_uattr_Njcomp();
double pairPot_integral_switch(double r,int icomp,int jcomp,int typePairPot);
extern int Type_pairPot;
double pairPot_ATT_noCS_switch(double r,int icomp,int jcomp,int typePairPot);
#define PI    M_PI
#define NCOMP_MAX 5
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
double StenTheta_uattr_sten_vol(int i,int j);
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
double StenTheta_uattr_sten_rad(int icomp,int jcomp);
