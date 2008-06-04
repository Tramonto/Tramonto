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
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
void sort_by_real(int nconv,int ncv,int ldv,double *d,double *v);
void sort_by_real(int nconv,int ncv,int ldv,double *d,double *v);
#if defined(LOCA_MF) && !defined(_CON_CONST_H_)
#include <mf.h>
#endif
void mass_matvec_mult_conwrap(double *x,double *y);
#if !defined(_CON_CONST_H_)
void mass_matvec_mult_conwrap(double *x,double *y);
#endif
void mass_matvec_mult_conwrap(double *x,double *y);
void mass_matvec_mult_conwrap(double *x,double *y);
void matvec_mult_conwrap(double *x,double *y);
#if !defined(_CON_CONST_H_)
void matvec_mult_conwrap(double *x,double *y);
#endif
void matvec_mult_conwrap(double *x,double *y);
void matvec_mult_conwrap(double *x,double *y);
double null_vector_resid(double r_val,double i_val,double *r_vec,double *i_vec,int mm_flag);
double null_vector_resid(double r_val,double i_val,double *r_vec,double *i_vec,int mm_flag);
double ltransnorm(double *x,double *scale_vec);
double ltransnorm(double *x,double *y);
double ip(double *x,double *y);
double ip(double *x,double *y);
double scaled_dp(double *x,double *y);
double scaled_dp(double *x,double *y);
double dp(double *x,double *y);
double dp(double *x,double *y);
void vec_copy(double *dx,double *dy);
void vec_copy(double *dx,double *dy);
void free_vec(double **ptr);
void free_vec(double **ptr);
double *alloc_vec();
double *alloc_vec();
void vec_init(double *u);
void vec_init(double *u);
void init_scale_utils(double *scale_vec);
void init_scale_utils(double *scale_vec);
void initialize_util_routines(int n_o,int n_t);
void initialize_util_routines(int n_o,int n_t);
extern double *SclVec;
extern double N_g;
extern int N_t;
extern int N_o;
double gsum_double_conwrap(double sum);
double gsum_double_conwrap(double sum);
#if !defined(_CON_CONST_H_)
double gsum_double_conwrap(double sum);
#endif
double gsum_double_conwrap(double sum);
double gsum_double_conwrap(double sum);
