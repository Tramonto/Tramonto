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
double StenDelta_Bond_GetWeightFromSten(double rsq,double R);
int StenDelta_Bond_NquadPtsGauss(double r);
extern int Ndim;
int StenDelta_Bond_NquadPtsBoundary();
extern int Ncomp;
int StenDelta_Bond_Njcomp();
double StenDelta_Bond_sten_vol(int icomp,int jcomp);
#define NCOMP_MAX 5
extern double Bond_ff[NCOMP_MAX][NCOMP_MAX];
double StenDelta_Bond_sten_rad(int icomp,int jcomp);
