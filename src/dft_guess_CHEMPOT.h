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
#define INIT_GUESS_FLAG  2
double load_nonlinear_transport_eqn(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
void node_box_to_ijk_box(int node_box,int *ijk_box);
#define NEQ_TYPE       12 
extern int Phys2Unk_last[NEQ_TYPE];
extern void *LinProbMgr_manager;
void calc_init_chem_pot(double **xInBox,double **xOwned);
#define NMER_MAX     200
extern double Betamu_chain_RTF[NMER_MAX];
extern double Betamu_chain_LBB[NMER_MAX];
#define NCOMP_MAX 5
extern double Betamu_RTF[NCOMP_MAX];
extern double Betamu_LBB[NCOMP_MAX];
#define POISSON        1
extern double Charge_f[NCOMP_MAX];
#define DENSITY        0
extern int Ipot_ff_c;
extern double VEXT_MAX;
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int *L2B_node;
extern int **Zero_density_TF;
extern int Unk2Comp[NMER_MAX];
#define DIFFUSION      6
extern int Phys2Unk_first[NEQ_TYPE];
extern int Nmissing_densities;
#define RESTART_FEWERCOMP  4
extern int Restart;
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
void node_to_ijk(int node,int *ijk);
extern int *L2G_node;
extern int Nnodes_per_proc;
extern double X_const_mu;
extern int Grad_dim;
extern double Size_x[NDIM_MAX];
extern int Npol_comp;
extern int Nseg_tot;
extern int Lseg_densities;
extern int Ncomp;
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_poly;
void setup_chem_pot(double **xOwned);
