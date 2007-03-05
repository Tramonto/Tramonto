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
#define PAIR_COULOMB       2
extern int Type_pairPot;
#define NONE       -1
#define NONE      -1
#define NONE -1
#define NONE        -1
extern int Type_attr;
#define NCOMP_MAX 5
extern double Charge_f[NCOMP_MAX];
double integrate_potential(double param1,double param2,double param3,int ngp,int ngpu,double *gp,double *gpu,double *gw,double *gwu,double *node_pos,double *node_pos_f);
extern int **Zero_density_TF;
extern int *L2B_node;
void find_images_coulomb(int idim,int *image,double **image_pos,double *node_image);
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
void node_to_position(int inode,double *NodePos);
int element_to_node(int ielement);
void set_gauss_quad(int ngp,double *gp,double *gw);
extern int Ndim;
int node_to_elem(int inode_all,int local_node,int *reflect_flag);
void safe_free(void **ptr);
void safe_free(void **ptr);
extern int *Comm_offset_node;
extern int *Comm_node_proc;
extern int *L2G_node;
extern double *Charge_vol_els;
extern int *B2L_node;
int element_box_to_node_box(int iel_box);
extern int Nelements_box;
extern int Nnodes;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Nelements;
extern int Ncomp;
extern int Nnodes_per_proc;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
extern double **Vext_coul;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
void setup_vext_coulomb_vol();
