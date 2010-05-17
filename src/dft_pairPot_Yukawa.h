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
double uYUKAWA_Integral(double r,int i,int j);
double uYUKAWA_ATT_CS(double r,int i,int j);
#define CORECONST_ZERO      1
double uYUKAWA_ATT_noCS(double r,int i,int j);
#define CORECONST_UCONST    0
extern int Type_CoreATT_CONST;
void uYUKAWA_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
double uYUKAWA_DERIV1D(double r,double x,double sigma,double eps,double rcut,double yukawaK);
#define NWALL_MAX_TYPE 50 
extern double YukawaK_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Cut_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
#define WALL_WALL   2
#define NCOMP_MAX 5
extern double YukawaK_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Cut_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
#define NWALL_MAX 600 
extern int WallType[NWALL_MAX];
extern double Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];
#define WALL_FLUID  1
extern double YukawaK_ff[NCOMP_MAX][NCOMP_MAX];
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
#define FLUID_FLUID 0
void uYUKAWA_CS_setparams(int context,int i,int j,double *param1,double *param2,double *param3,double *param4);
double uYUKAWA_CS(double r,double sigma,double eps,double rcut,double yukawaK);
