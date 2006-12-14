/* This file was automatically generated.  Do not edit! */
int node_box_to_elem_box_reflect(int inode_box,int local_node,int *reflect_flag);
int node_box_to_elem_box(int inode_box,int local_node);
int el_to_el_box(int iel);
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
extern int Elements_plane_box;
#define NDIM_MAX  3
extern int Elements_x_box[NDIM_MAX];
int element_box_to_node_box(int iel_box);
int el_box_to_el(int iel_box);
int unk_to_unk_box(int i);
extern int Nunk_per_node;
int unk_box_to_unk(int i_box);
int node_to_node_box_no_bound(int inode);
int node_to_node_box(int inode);
int node_box_to_node(int inode_box);
int ijk_box_to_node_box(int *ijk_box);
void node_box_to_ijk_box(int node_box,int *ijk_box);
int position_to_node(double *NodePos);
extern double Esize_x[NDIM_MAX];
extern double Size_x[NDIM_MAX];
void node_to_position(int inode,double *NodePos);
int node_to_elem_v2(int inode_all,int local_node);
int node_to_elem_return_dim(int inode_all,int local_node,int *reflect_flag,int *idim_return,int *iside,int *periodic_flag);
int node_to_elem(int inode_all,int local_node,int *reflect_flag);
void node_to_ijk(int node,int *ijk);
extern int Nnodes_per_el_V;
void element_to_nodes(int ielement,int *nodes);
extern int Nodes_plane;
int map_0th_plane(int i,int Nplane);
extern int Elements_plane;
extern int Elements_x[NDIM_MAX];
int element_to_node(int ielement);
extern int Nodes_plane_box;
extern int Min_IJK[3];
extern int Pflag[3];
extern int Max_IJK[3];
extern int Max_IJK_box[3];
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
void ijk_box_to_ijk(int *ijk_box,int *ijk);
void ijk_to_ijk_box_no_bound(int *ijk,int *ijk_box);
extern int Nodes_x_box[NDIM_MAX];
extern int Min_IJK_box[3];
void ijk_to_ijk_box(int *ijk,int *ijk_box);
int ijk_to_node(int *ijk);
#define LAST_NODE    3
#define REFLECT      2
extern int Nodes_x[NDIM_MAX];
#define PERIODIC     1
#define IN_WALL     -1
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Lsteady_state;
extern int Iliq_vap;
extern int Nwall;
#define IN_BULK      0
extern int Type_bc[NDIM_MAX][2];
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Ndim;
int offset_to_node(int *inode_ijk,int *offset_node,int *reflect_flag);
int round_to_int(double x);
