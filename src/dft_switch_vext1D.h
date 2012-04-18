/* This file was automatically generated.  Do not edit! */
double Vextderiv_LINEAR(double x,int icomp,int iwall_type);
double Vextderiv_LJ7YukawaSum(double r,int icomp,int iwall_type);
double Vextderiv_EXP_ATT(double x,int icomp,int iwall_type);
double Vextderiv_REPULSIVE9(double x,int icomp,int iwall_type);
double Vextderiv_LJ9_3_v2(double x,int icomp,int iwall_type);
double Vextderiv_LJ9_3(double x,int icomp,int iwall_type);
double Vext_1D_dash(double x,int icomp,int iwall_type);
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
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
double Vext_LINEAR_noCS(double x,int icomp,int iwall_type);
#define LINEAR_noCS       6 
double Vext_LJ7YukawaSum(double r,int icomp,int iwall_type);
#define R7_YUKAWA_SUM_CS  7 
double Vext_EXP_ATT_noCS(double x,int icomp,int iwall_type);
#define EXP_ATT_noCS      5
double Vext_REPULSIVE9_noCS(double x,int icomp,int iwall_type);
#define REPULSIVE9_noCS   4
double Vext_LJ9_3_shiftX_CS(double x,int icomp,int iwall_type);
#define LJ9_3_shiftX_CS   3
double Vext_LJ9_3_noCS(double x,int icomp,int iwall_type);
#define LJ9_3_noCS        2
double Vext_LJ9_3_v2_CS(double x,int icomp,int iwall_type);
#define LJ9_3_v2_CS       1
double Vext_LJ9_3_CS(double x,int icomp,int iwall_type);
#define LJ9_3_CS          0
#define NWALL_MAX_TYPE 20 
extern int Vext_PotentialID[NWALL_MAX_TYPE];
double Vext_1D(double x,int icomp,int iwall_type);
