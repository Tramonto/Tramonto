/* This file was automatically generated.  Do not edit! */
/*void print_to_file_comp(FILE *fp,int icomp,double val,char *var_label,int first);
void print_to_file(FILE *fp,double val,char *var_label,int first);*/
void print_to_screen_comp(int icomp,double val,char *var_label);
void print_to_screen(double val,char *var_label);
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
#define LAST_NODE    3
#define IN_BULK      0
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
int el_to_el_box(int iel);
int node_to_elem(int inode_all,int local_node,int *reflect_flag);
extern int Imax;
extern int List[2];
extern int **Nel_hit2;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
extern int **Nel_hit;
void setup_integrals();
#define NWALL_MAX 600 
extern int Link[NWALL_MAX];
extern double **S_area_tot;
extern int Nlink;
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int Lcount_reflect;
#define REFLECT      2
void setup_domain_multipliers();
#define PERIODIC     1
extern int Type_bc[NDIM_MAX][2];
extern int Nodes_x[NDIM_MAX];
extern int *L2G_node;
void node_to_ijk(int node,int *ijk);
void integrateOverSurface(double(*fp_integrand)(int,int,int,double **),int iunk,double **x,double *profile);
void integrateInSpace_SumInComp(double(*fp_integrand)(int,int,double **),int **nelhit,double **x,double *profile);
extern double Fac_area;
extern double Fac_vol;
extern double Area;
extern int Lper_area;
double gsum_double(double c);
extern double Vol_el;
extern int *L2B_node;
void integrateInSpace(double(*fp_integrand)(int,int,double **),int iunk,int **nelhit,double **x,double *profile);
extern int Nnodes;
#define GLOBAL 2
#define BOX 0
extern int Nnodes_per_proc;
#define LOCAL_N 1
extern int Nunk_per_node;
#define MATRIX_FILL_NODAL 1   /* set to zero for physics based ordering */
int loc_find(int iunk,int inode,int flag);
#define CMS_G          2 
#define NMER_MAX     100
extern double BondWTC_RTF[NMER_MAX *NMER_MAX];
extern double BondWTC_LBB[NMER_MAX *NMER_MAX];
extern double BondWTC_b[NMER_MAX *NMER_MAX];
#define BONDWTC       7
extern double Xi_cav_RTF[4];
extern double Xi_cav_LBB[4];
extern double Xi_cav_b[4];
#define CAVWTC     6
#define NCOMP_MAX 5
extern double Betamu_RTF[NCOMP_MAX];
extern double Betamu_LBB[NCOMP_MAX];
extern double VEXT_MAX;
#define DIFFUSION      5
extern double Elec_pot_RTF;
extern double Elec_pot_LBB;
#define POISSON        3
extern double Rhobar_b_RTF[10];
extern double Rhobar_b_LBB[10];
extern double Rhobar_b[10];
#define HSRHOBAR       4
extern double Rho_b_RTF[NCOMP_MAX];
extern double Rho_coex[2];
extern double Rho_b_LBB[NCOMP_MAX];
extern double Rho_seg_RTF[NMER_MAX];
extern double Rho_seg_LBB[NMER_MAX];
extern int Lsteady_state;
extern double Rho_seg_b[NMER_MAX];
extern int Unk2Phys[3 *NCOMP_MAX+NMER_MAX+NMER_MAX *NMER_MAX+13];
extern int *Mesh_coarsen_flag;
extern int Nwall;
extern int Mesh_coarsening;
extern int Nzone;
extern int Coarser_jac;
extern int WallType[NWALL_MAX];
extern int **Lsemiperm;
extern int **Wall_elems;
int node_box_to_elem_box_reflect(int inode_box,int local_node,int *reflect_flag);
extern int Nnodes_per_el_V;
double resid_and_Jac_sten_fill(int sten_type,double **x,int iunk,int junk,int icomp,int jcomp,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag,int jzone_flag,double(*fp_prefactor)(int,int,int *),double(*fp_resid)(int,int,double **),double(*fp_jacobian)(int,int,double **));
extern void *LinProbMgr_manager;
extern double Rho_b[NCOMP_MAX];
#define CMS_FIELD      1
#define POLYMER_CR     4
#define THETA_CHARGE   3
#define U_ATTRACT      2
extern int *Pol_Sym_Seg;
extern int Nseg_tot;
#define FLAG_PBELEC -777
#define FLAG_BULK   -888
int find_jzone(int izone,int inode_box);
double resid_and_Jac_sten_fill_sum_Ncomp(int sten_type,double **x,int iunk,int loc_inode,int inode_box,int izone,int *ijk_box,int resid_only_flag,int jzone_flag,double(*fp_prefactor)(int,int,int *),double(*fp_resid)(int,int,double **),double(*fp_jacobian)(int,int,double **));
extern double Temporary_sum;
double constant_boundary(int iunk,int jnode_box);
extern int Nnodes_box;
double HW_boundary_weight(int icomp,int ilist,double *hw_weight,int inode_box,int *reflect_flag);
extern int **Nodes_2_boundary_wall;
extern int Lhard_surf;
extern int **Zero_density_TF;
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
extern int Nlists_HW;
#define NEQ_TYPE       8
extern int Phys2Unk_last[NEQ_TYPE];
#define DENSITY        0
extern int Phys2Unk_first[NEQ_TYPE];
extern int Unk2Comp[NMER_MAX];
#define WTC          3
extern int Type_poly;
void node_box_to_ijk_box(int node_box,int *ijk_box);
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
void int_stencil(double **x,int inode_box,int iunk,int sten_type);
extern double Esize_x[NDIM_MAX];
extern int Ndim;
extern int Ncomp;
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
struct Stencil_Struct {
  int        Length;      /* Number of nodes that interact with current 
                             node through this stencil                    */
  int      **Offset;      /* 2D array to be set to size [Length][Ndim] that
                             gives relative position of interacting node  */
  double    *Weight;      /* 1D array of size [Length] that gives weight
                             of interacting node to total stencil         */
  double   **HW_Weight;   /* 2D array of size [Length][Nnodes_per_el] that
                             holds the weights based on what element they
                             are being contributed from. Only used for Hard
                             Walls when stencil point is a boundary node  */
};
double int_stencil_bulk(int sten_type,int icomp,int jcomp,double(*fp_integrand)(double,int,int));
