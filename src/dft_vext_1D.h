/* This file was automatically generated.  Do not edit! */
double Vextderiv_LINEAR(double x,int icomp,int iwall_type);
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
#define NWALL_MAX 600 
extern double Elec_param_w[NWALL_MAX];
#define NCOMP_MAX 5
extern double Charge_f[NCOMP_MAX];
extern double Temp_elec;
double Vext_LINEAR_noCS(double x,int icomp,int iwall_type);
double Vextderiv_EXP_ATT(double x,int icomp,int iwall_type);
double Vext_EXP_ATT_noCS(double x,int icomp,int iwall_type);
double Vextderiv_REPULSIVE9(double x,int icomp,int iwall_type);
double Vext_REPULSIVE9_noCS(double x,int icomp,int iwall_type);
extern double VEXT_MAX;
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
double Vext_LJ9_3_shiftX_CS(double x,int icomp,int iwall_type);
double Vext_LJ9_3_noCS(double x,int icomp,int iwall_type);
double Vextderiv_LJ9_3_v2(double x,int icomp,int iwall_type);
double Vext_LJ9_3_v2_CS(double x,int icomp,int iwall_type);
double Vextderiv_LJ9_3(double x,int icomp,int iwall_type);
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
#define PI    3.141592653589793238462643383279502884197169399375
#define NWALL_MAX_TYPE 50 
extern double Cut_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Rho_w[NWALL_MAX_TYPE];
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
double Vext_LJ9_3_CS(double x,int icomp,int iwall_type);
