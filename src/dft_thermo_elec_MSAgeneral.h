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
#define VERBOSE      3 
extern int Iwrite;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void precalc_GENmsa_params(double *rho,double *x_msa,double *n_msa,double gamma);
#define PI    3.141592653589793238462643383279502884197169399375
double deltaC_GENERAL_MSA_int(double r,int i,int j);
#define NCOMP_MAX 5
extern double MSAgen_term4[NCOMP_MAX][NCOMP_MAX];
extern double MSAgen_term3[NCOMP_MAX][NCOMP_MAX];
extern double MSAgen_term2[NCOMP_MAX][NCOMP_MAX];
extern double MSAgen_term1[NCOMP_MAX][NCOMP_MAX];
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
extern double Charge_f[NCOMP_MAX];
extern double Temp_elec;
extern double X_MSA[NCOMP_MAX];
extern double Gamma_MSA;
extern double N_MSA[NCOMP_MAX];
extern double HS_diam[NCOMP_MAX];
double deltaC_GENERAL_MSA(double r,int i,int j);
#define THETA_CR_GENERAL_MSA  7
double int_stencil_bulk(int sten_type,int icomp,int jcomp,double(*fp_integrand)(double,int,int));
extern int Ncomp;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern double *Deltac_b;
void chempot_ELEC_MSA_GENERAL(double *rho);
