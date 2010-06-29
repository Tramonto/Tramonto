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
int StenTheta_uattr_NquadPtsBoundary();
extern int Ncomp;
int StenTheta_uattr_Njcomp();
#define PAIR_SW		      5
double pairPot_integral_switch(double r,int icomp,int jcomp,int typePairPot);
double pairPot_ATT_noCS_switch(double r,int icomp,int jcomp,int typePairPot);
#define PI    3.141592653589793238462643383279502884197169399375
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Type_pairPot;
void pairPot_InnerCore_switch(int icomp,int jcomp,int typePairPot,double *rCore_left,double *rCore_right,double *epsCore);
double StenTheta_uattr_sten_vol(int i,int j);
#define NCOMP_MAX 5
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
double StenTheta_uattr_sten_rad(int icomp,int jcomp);
