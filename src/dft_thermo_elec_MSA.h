/* This file was automatically generated.  Do not edit! */
double deltaC_MSA_int(double r,int i,int j);
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
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
#define NCOMP_MAX 5
extern double HS_diam[NCOMP_MAX];
extern double Charge_f[NCOMP_MAX];
extern double Rho_b[NCOMP_MAX];
extern double Temp_elec;
#define PI    3.141592653589793238462643383279502884197169399375
double deltaC_MSA(double r,int i,int j);
#define THETA_CR_RPM_MSA      3
double int_stencil_bulk(int sten_type,int icomp,int jcomp,double(*fp_integrand)(double,int,int));
extern int Ncomp;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern double *Deltac_b;
void chempot_ELEC_MSA_RPM(double *rho);
