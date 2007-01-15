/* This file was automatically generated.  Do not edit! */
double int_stencil_bulk(int sten_type,int icomp,int jcomp,double(*fp_integrand)(double,int,int));
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
extern int Ncomp;
#define NCOMP_MAX 5
extern double Rho_b[NCOMP_MAX];
double integrand_CMS_freen_bulk(int iunk,int inode_box,double **x);
extern int Nmer[NCOMP_MAX];
#define POLYMER_CR     4
double int_stencil(double **x,int inode_box,int iunk,int sten_type);
#define NBLOCK_MAX   5
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
#define DENSITY        0
#define NEQ_TYPE       8
extern int Phys2Unk_first[NEQ_TYPE];
double integrand_CMS_freen(int iunk,int inode_box,double **x);
