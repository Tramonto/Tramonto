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
extern double **Vext;
void setup_exp_density_with_profile(double **xInBox);
#define NCOMP_MAX 5
#define NSTEPS_MAX 10
extern double Rho_step[NCOMP_MAX][NSTEPS_MAX];
#define CHOP_RHO_STEP    4
#define CHOP_RHO         3
extern double Xstart_step[NSTEPS_MAX];
extern double **X_wall;
extern int Nwall;
void chop_profile(double **xInBox,int guess_type);
void check_zero_densities_owned(double **xOwned);
extern double Rho_b[NCOMP_MAX];
extern double VEXT_MAX;
#define NMER_MAX     100
extern double Rho_seg_b[NMER_MAX];
extern int **Zero_density_TF;
extern void *LinProbMgr_manager;
extern int *L2B_node;
extern int Nnodes_per_proc;
void communicate_to_fill_in_box_values(double **xInBox);
extern int *B2G_node;
extern int Nnodes_box;
int locate_inode_old(int *ijk);
extern int Pos_new_nodes;
void node_to_ijk(int node,int *ijk);
#define NWALL_MAX_TYPE 50 
extern double Del_1[NWALL_MAX_TYPE];
extern int Plane_new_nodes;
#define NDIM_MAX  3
extern int Nodes_x_old[NDIM_MAX];
extern int Nodes_x[NDIM_MAX];
#define NEQ_TYPE       13 
extern int Phys2Unk_last[NEQ_TYPE];
extern int Nrho_bar_s;
extern int Nmer_comp[NCOMP_MAX];
extern int Phys2Unk_first[NEQ_TYPE];
extern double Temp;
extern double Mass[NCOMP_MAX];
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern int LDeBroglie;
extern double Scale_fac_WJDC[NCOMP_MAX][NCOMP_MAX];
extern int **Nseg_type_pol;
extern int Npol_comp;
extern int Unk2Comp[NMER_MAX];
int ijk_to_node(int *ijk);
extern double Esize_x[NDIM_MAX];
int round_to_int(double x);
#define G_CHAIN       11 
#define CALC_RHOBAR_ONLY  2
#define BULK              0
#define CALC_RHOBAR_AND_G 3
#define CALC_ALL_FIELDS   1
extern int Iguess_fields;
#define VERBOSE      3 
#define MF_VARIABLE  2
extern int Type_attr;
extern int L_HSperturbation;
#define CMS_SCFT     1
#define CMS          0
#define NONE       -1
#define NONE      -1
#define NONE        -1
#define NONE        -1
extern int Type_coul;
#define DIFFUSIVE_INTERFACE 1
extern int Type_interface;
#define NO_SCREEN    2 
extern int Iwrite;
#define DIFFUSION      6
extern int Nbonds;
#define BONDWTC        5
#define CAVWTC         4
#define WJDC3        5 
#define WJDC2        4 
#define WJDC         3
extern int Type_poly;
#define WJDC_FIELD     8
#define SCF_CONSTR	   9
#define SCF_FIELD	  10
#define CMS_FIELD      7
extern int Nrho_bar;
extern int Ndim;
#define HSRHOBAR       2
#define POISSON        1
#define MF_EQ          3
extern int Nseg_tot;
extern int Lseg_densities;
extern int Nmissing_densities;
extern int Ncomp;
#define RESTART_FEWERCOMP  4
void check_zero_densities(double **xInBox);
void communicate_profile(double *x_new,double **xInBox);
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define DENSITY        0
void safe_free(void **ptr);
void safe_free(void **ptr);
extern int Iguess;
extern int Nunknowns;
void shift_the_profile(double *x_new,double fac,double *xold);
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Restart_field[NEQ_TYPE];
void read_in_a_file(int guess_type,char *filename);
extern double *X_old;
extern double *X2_old;
#define RESTART_1DTOND     5
extern int Restart;
int find_length_of_file(char *filename);
extern int Nodes_old;
#define BINODAL_FLAG  -1325  /* to let initial guess routine know we need to fill X2 */
extern int Lbinodal;
extern int Imain_loop;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Nunk_per_node;
extern int Nnodes;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
void guess_restart_from_files(int start_no_info,int guess_type,double **xInBox);
