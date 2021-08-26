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
extern double Rho_b[NCOMP_MAX];
#define NMER_MAX     200
extern double Rho_seg_b[NMER_MAX];
extern int Lseg_densities;
#define DENSITY        0
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Icomp_to_polID[NCOMP_MAX];
extern int Grafted[NCOMP_MAX];
#define WJDC3        5 
extern int Type_poly;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Grafted_Logical;
double integrand_ideal_gas_freen_bulk(int iunk,int inode_box,double **x);
#define DENSITY_MIN  1.e-20
double integrand_ideal_gas_freen(int iunk,int inode_box,double **x);
