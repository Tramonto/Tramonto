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
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define NCOMP_MAX 5
extern double Charge_f[NCOMP_MAX];
extern int Ipot_ff_c;
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
extern double D_coef[NCOMP_MAX];
extern int Linear_transport;
extern double Velocity;
int ijk_to_node(int *ijk);
extern int Nunk_per_node;
#define DENSITY        0
#define DIFFUSION      5
#define NEQ_TYPE       8
extern int Phys2Unk_first[NEQ_TYPE];
extern int Ncomp;
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
extern int Nodes_x[NDIM_MAX];
extern int Ndim;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
void node_to_ijk(int node,int *ijk);
extern int Nnodes;
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern double *X_old;
void calc_flux(FILE *fp,char *output_flux,double *X_old);
