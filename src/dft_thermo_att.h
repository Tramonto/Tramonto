/* This file was automatically generated.  Do not edit! */
double dmu_drho_att(double *rho);
double dp_drho_att(double *rho);
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
#define NCOMP_MAX 5
extern double Betamu_att[NCOMP_MAX];
void chempot_att(double *rho);
double pressure_att(double *rho);
#define THETA_PAIRPOT_RCUT    2
double int_stencil_bulk(int sten_type,int icomp,int jcomp,double(*fp_integrand)(double,int,int));
extern int Icomp_to_polID[NCOMP_MAX];
extern int Grafted[NCOMP_MAX];
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Grafted_Logical;
extern double Avdw[NCOMP_MAX][NCOMP_MAX];
extern int Ncomp;
void calc_Avdw_att();
void ATT_thermo_precalc();
