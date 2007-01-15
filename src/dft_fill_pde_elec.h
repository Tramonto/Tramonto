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
extern double Area_surf_el[3];
extern double **Charge_w_sum_els;
#define CONST_CHARGE     2
extern double *Dielec;
#define NCOMP_MAX 5
extern double Pol[NCOMP_MAX];
extern int Lpolarize[NCOMP_MAX];
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
extern double Charge_f[NCOMP_MAX];
extern int **Lsemiperm;
extern int **Wall_elems;
#define NMER_MAX     100
extern int Unk2Comp[NMER_MAX];
extern int Lseg_densities;
#define NEQ_TYPE       8
extern int Phys2Unk_last[NEQ_TYPE];
#define DENSITY        0
#define KAPPA_H2O 78.5
extern int **Zero_density_TF;
extern int Ncomp;
#define POISSON        3
extern int Phys2Unk_first[NEQ_TYPE];
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
#if defined(DEC_ALPHA)
#define POW_INT powii
#endif
#if !(defined(DEC_ALPHA))
#define POW_INT (int)pow
#endif
extern double Temp_elec;
#define PI    M_PI
extern double *Charge_vol_els;
extern int Vol_charge_flag;
int el_to_el_box(int iel);
int node_to_elem(int inode_all,int local_node,int *reflect_flag);
extern int Nnodes_per_el_V;
#define NWALL_MAX 600 
extern double Elec_param_w[NWALL_MAX];
#define CONST_POTENTIAL  1
extern int WallType[NWALL_MAX];
#define NWALL_MAX_TYPE 50 
extern int Type_bc_elec[NWALL_MAX_TYPE];
extern int Nlists_HW;
extern int **Nodes_2_boundary_wall;
void set_fem_1el_weights(double **wt_lp_1el_ptr,double **wt_s_1el_ptr,int ***elem_permute);
extern int Ndim;
double load_poisson_bc(int iunk,int loc_inode,int inode_box);
double load_poissons_eqn(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
double load_polarize_poissons_eqn(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
#define POLARIZE   2
extern int Type_coul;
extern double Elec_pot_RTF;
extern void *LinProbMgr_manager;
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
extern int Lsteady_state;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
double load_poisson_control(int iunk,int loc_inode,int inode_box,int *ijk_box,double **x,int resid_only_flag);
