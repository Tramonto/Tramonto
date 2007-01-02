/* This file was automatically generated.  Do not edit! */
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
#define PI    M_PI
double uCOULOMB_Integral(double r,int i,int j);
#define NCOMP_MAX 5
extern double Charge_f[NCOMP_MAX];
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
double uCOULOMB_ATT_noCS(double r,int i,int j);
double uCOULOMB_DERIV1D(double r,double x,double z1,double z2);
extern double Temp_elec;
double uCOULOMB(double r,double z1,double z2);
