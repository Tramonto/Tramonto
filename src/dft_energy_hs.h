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
extern int Ndim;
extern double Rhobar_b[10];
extern double Rhobar_b_RTF[10];
extern int Lsteady_state;
extern int Nrho_bar_s;
double integrand_hs_freen_bulk(int iunk,int inode_box,double **x);
double phispt_switch(double *n);
void solutionVec_to_nOrdering(double *rhoBar_SVOrdering,double *n);
#define NEQ_TYPE       10 
extern int Phys2Unk_first[NEQ_TYPE];
#define HSRHOBAR       4
extern int Phys2Nunk[NEQ_TYPE];
#define NDIM_MAX  3
double integrand_hs_freen(int iunk,int inode_box,double **x);
