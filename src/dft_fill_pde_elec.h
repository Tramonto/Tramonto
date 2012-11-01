/* This file was automatically generated.  Do not edit! */
double load_poisson_bc(int iunk,int loc_inode,int inode_box,int resid_only_flag);
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
extern double *Dielec;
void ijk_box_to_ijk(int *ijk_box,int *ijk);
extern double Area_surf_el[3];
extern double **Charge_w_sum_els;
#define CONST_CHARGE     2
#define NCOMP_MAX 5
extern double Pol[NCOMP_MAX];
extern int Lpolarize[NCOMP_MAX];
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
#define KAPPA_H2O 78.5
extern double Charge_f[NCOMP_MAX];
extern int **Lsemiperm;
extern int **Wall_elems;
#define NMER_MAX     200
extern int Unk2Comp[NMER_MAX];
extern int Lseg_densities;
#define NEQ_TYPE       12 
extern int Phys2Unk_last[NEQ_TYPE];
#define DENSITY        0
extern int **Zero_density_TF;
extern int Ncomp;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
#if defined(DEC_ALPHA)
#define POW_INT powii
#endif
#if !(defined(DEC_ALPHA))
#define POW_INT (int)pow
#endif
#define POISSON        1
extern int Phys2Unk_first[NEQ_TYPE];
extern double Temp_elec;
#define PI    3.141592653589793238462643383279502884197169399375
extern double *Charge_vol_els;
extern int Vol_charge_flag;
int el_to_el_box(int iel);
int node_to_elem(int inode_all,int local_node,int *reflect_flag);
extern int Nnodes_per_el_V;
#define NWALL_MAX 600 
extern double Elec_param_w[NWALL_MAX];
#define CONST_POTENTIAL  1
extern int WallType[NWALL_MAX];
#define NWALL_MAX_TYPE 20 
extern int Type_bc_elec[NWALL_MAX_TYPE];
extern int Nlists_HW;
extern int **Nodes_2_boundary_wall;
void set_fem_1el_weights(double **wt_lp_1el_ptr,double **wt_s_1el_ptr,int ***elem_permute);
extern int Ndim;
double load_poissons_eqn(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
double load_polarize_poissons_eqn(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
#define POLARIZE       3
extern int Type_coul;
extern double Elec_pot_RTF;
extern int *B2G_node;
extern int Nnodes;
extern int Solver_Unk[3 *NCOMP_MAX+2 *NMER_MAX+NMER_MAX *NMER_MAX+13];
extern int *L2G_node;
extern double **Array_test;
#define FILES_DEBUG_MATRIX 3 
extern int Iwrite_files;
extern void *LinProbMgr_manager;
#define CALC_RESID_ONLY  3
#define INIT_GUESS_FLAG  2
extern double Elec_pot_LBB;
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern double X_const_mu;
extern double Size_x[NDIM_MAX];
extern int Grad_dim;
int node_box_to_node(int inode_box);
void node_to_position(int inode,double *NodePos);
#define UNIFORM_INTERFACE  0
extern int Type_interface;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
double load_poisson_control(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
