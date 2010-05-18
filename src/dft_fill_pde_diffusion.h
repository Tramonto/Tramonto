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
extern int **Nodes_2_boundary_wall;
void set_fem_1el_weights(double **wt_lp_1el_ptr,double **wt_s_1el_ptr,int ***elem_permute);
double load_linear_transport_eqn(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
extern double Velocity;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
#if defined(DEC_ALPHA)
#define POW_INT powii
#endif
#if !(defined(DEC_ALPHA))
#define POW_INT (int)pow
#endif
#define DENSITY        0
extern double *Area_IC;
#define NWALL_MAX 600 
extern int WallType[NWALL_MAX];
extern int **Lsemiperm;
int el_to_el_box(int iel);
extern int **Wall_elems;
int node_to_elem_v2(int inode_all,int local_node);
extern int Nnodes_per_el_V;
#define NCOMP_MAX 5
extern double Betamu_RTF[NCOMP_MAX];
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
extern double Betamu_LBB[NCOMP_MAX];
extern double X_const_mu;
extern double Esize_x[NDIM_MAX];
extern int Grad_dim;
extern void *LinProbMgr_manager;
extern double VEXT_MAX;
extern int **Zero_density_TF;
int node_box_to_node(int inode_box);
void node_to_ijk(int node,int *ijk);
extern int Nlists_HW;
#define DIFFUSION      6
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern int *B2G_node;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
void basis_fn_calc(double **phi,double ***grad_phi,double *evol);
extern int Ndim;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
double load_nonlinear_transport_eqn(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
