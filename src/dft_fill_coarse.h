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
extern int Nzone;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
extern int *Mesh_coarsen_flag;
void node_box_to_ijk_box(int node_box,int *ijk_box);
double constant_boundary(int iunk,int jnode_box);
double load_coarse_variable(double **x,int jnode_box,double fac,int iunk,int loc_inode,int resid_only_flag);
void locate_neighbor_unks(double **x,int iunk,int loc_inode,int node_box,double fac,double *resid,int resid_only_flag);
#define POISSON        1
#define DENSITY        0
#define NCOMP_MAX 5
#define NMER_MAX     200
extern int Unk2Phys[3 *NCOMP_MAX+2 *NMER_MAX+NMER_MAX *NMER_MAX+13];
double load_coarse_node_Ndim(int loc_inode,int inode_box,int iunk,double **x,int resid_only_flag);
extern void *LinProbMgr_manager;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int *B2L_node;
int ijk_box_to_node_box(int *ijk_box);
extern int Grad_dim;
extern int Ndim;
void load_coarse_node_1dim(int loc_inode,int inode_box,int *ijk_box,int iunk,double **x,int resid_only_flag);
