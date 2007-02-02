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
extern double **Vext;
void setup_exp_density_with_profile(double **xOwned);
#define NCOMP_MAX 5
#define NSTEPS_MAX 10
extern double Rho_step[NCOMP_MAX][NSTEPS_MAX];
#define CHOP_RHO_STEP    7
#define CHOP_RHO         4
#define CHOP_RHO_V       6
extern double Rho_coex[2];
#define CHOP_RHO_L       5
extern double Xstart_step[NSTEPS_MAX];
extern double **X_wall;
extern int Nwall;
void chop_profile(double **xOwned,int iguess);
extern double VEXT_MAX;
#define NMER_MAX     100
extern double Rho_seg_b[NMER_MAX];
extern int **Zero_density_TF;
extern int Unk2Comp[NMER_MAX];
extern int *L2B_node;
extern int *L2G_node;
extern int Nnodes_per_proc;
extern double Rho_b[NCOMP_MAX];
extern double Rho_max;
extern int Unk2Phys[3 *NCOMP_MAX+NMER_MAX+NMER_MAX *NMER_MAX+13];
int locate_inode_old(int *ijk);
extern int Pos_new_nodes;
void node_to_ijk(int node,int *ijk);
#define NWALL_MAX_TYPE 50 
extern double Del_1[NWALL_MAX_TYPE];
extern int Plane_new_nodes;
#define NDIM_MAX  3
extern int Nodes_x_old[NDIM_MAX];
extern int Nodes_x[NDIM_MAX];
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern int *Unk_to_Seg;
extern int *Unk_to_Poly;
#define NEQ_TYPE       8
extern int Phys2Unk_last[NEQ_TYPE];
extern int Phys2Unk_first[NEQ_TYPE];
extern int Nrho_bar_s;
int ijk_to_node(int *ijk);
extern double Esize_x[NDIM_MAX];
int round_to_int(double x);
#define CMS_G          2 
extern int Ndim;
extern int L_HSperturbation;
#define CMS_SCFT     1
#define CMS          0
extern int Type_poly;
#define NONE       -1
#define NONE      -1
#define NONE -1
#define NONE        -1
extern int Type_coul;
extern int Lsteady_state;
#define NO_SCREEN    2 
extern int Iwrite;
#define DIFFUSION      5
extern int Nbonds;
#define BONDWTC       7
#define CAVWTC     6
#define CMS_FIELD      1
extern int Nrho_bar;
#define HSRHOBAR       4
#define POISSON        3
extern int Ncomp;
extern int Nseg_tot;
extern int Lseg_densities;
void check_zero_densities(double **xOwned);
void communicate_profile(double *x_new,double **xOwned);
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
extern int Iguess1;
extern int Nunknowns;
void shift_the_profile(double *x_new,double fac);
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Restart_field[NEQ_TYPE];
void read_in_a_file(int iguess,char *filename);
extern double *X_old;
extern double *X2_old;
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
void guess_restart_from_files(int start_no_info,int iguess,double **xOwned);
