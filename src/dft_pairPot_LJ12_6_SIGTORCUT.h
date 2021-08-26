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
#define NCOMP_MAX 6
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
double uLJ12_6_ATT_SIGTORCUT_CS(double r,int i,int j);
double uLJ12_6_ATT_SIGTORCUT_noCS(double r,int i,int j);
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
void uLJ12_6_SIGTORCUT_InnerCore(int i,int j,double *rCore_left,double *rCore_right,double *epsCore);
