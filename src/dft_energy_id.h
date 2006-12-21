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
#define NCOMP_MAX 5
extern double Rho_b[NCOMP_MAX];
#define DENSITY        0
#define NEQ_TYPE       8
extern int Phys2Unk_first[NEQ_TYPE];
#define NMER_MAX     100
extern double Rho_seg_b[NMER_MAX];
#define WTC          3
extern int Type_poly;
double integrand_ideal_gas_freen_bulk(int iunk,int inode_box,double **x);
#define DENSITY_MIN  1.e-20
double integrand_ideal_gas_freen(int iunk,int inode_box,double **x);
