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
#define FMT3       2
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
#define FMT2       1
#define PI    M_PI
#define FMT1       0
extern int Type_func;
extern int Nrho_bar_s;
extern int Ndim;
#define NDIM_MAX  3
extern double Rhobar_b[10];
extern double Rhobar_b_RTF[10];
extern int Lsteady_state;
double integrand_hs_freen_bulk(int iunk,int inode_box,double **x);
double phispt(double *rho_bar);
#define NEQ_TYPE       8
extern int Phys2Unk_first[NEQ_TYPE];
#define HSRHOBAR       4
extern int Phys2Nunk[NEQ_TYPE];
double integrand_hs_freen(int iunk,int inode_box,double **x);
