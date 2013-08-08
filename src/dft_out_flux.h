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
#define NCOMP_MAX 5
extern double D_coef[NCOMP_MAX];
extern int Linear_transport;
extern double Velocity;
int ijk_to_node(int *ijk);
extern int Nunk_per_node;
#define DENSITY        0
#define DIFFUSION      6
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Ncomp;
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Nodes_x[NDIM_MAX];
extern int Ndim;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
void node_to_ijk(int node,int *ijk);
extern int Nnodes;
#define FILENAME_LENGTH 4096
extern char Outpath_array[FILENAME_LENGTH];
extern double *X_old;
void calc_flux(FILE *fp,char *output_flux,double *X_old);
