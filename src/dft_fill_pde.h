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
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void basis_fn_calc(double **phi,double ***grad_phi,double *evol);
extern int Nnodes_per_el_V;
void set_fem_1el_weights(double **wt_lp_1el_ptr,double **wt_s_1el_ptr,int ***elem_permute);
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern int Ndim;
#if defined(DEC_ALPHA)
#define POW_INT powii
#endif
#if !(defined(DEC_ALPHA))
#define POW_INT (int)pow
#endif
void set_fem_weights(double **wt_laplace_ptr,double **wt_source_ptr);
