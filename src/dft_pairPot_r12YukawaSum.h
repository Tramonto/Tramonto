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
double ur12andYUKAWA_Integral(double r,int i,int j);
double ur12andYUKAWA_ATT_CS(double r,int i,int j);
#define CORECONST_ZERO      1
double ur12andYUKAWA_ATT_noCS(double r,int i,int j);
#define CORECONST_UCONST    0
extern int Type_CoreATT_CONST;
#define ATTCORE_SIGTOUMIN   3
#define NCOMP_MAX 5
extern double Rzero_ff[NCOMP_MAX][NCOMP_MAX];
#define ATTCORE_UCSZERO     2
extern double Rmin_ff[NCOMP_MAX][NCOMP_MAX];
#define ATTCORE_UMIN        1
#define ATTCORE_SIGMA       0
extern int Type_CoreATT_R;
void ur12andYUKAWA_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
double ur12andYUKAWA_DERIV1D(double r,double x,double sigma,double eps,double rcut,double yukawaK);
#define NWALL_MAX_TYPE 50 
extern double YukawaK_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Cut_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
#define WALL_WALL   2
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
void ur12andYUKAWA_CS_setparams(int context,int i,int j,double *param1,double *param2,double *param3,double *param4);
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
double ur12andYUKAWA_CS(double r,double sigma,double eps,double rcut,double yukawaK);
