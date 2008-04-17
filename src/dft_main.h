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
extern double **Charge_x;
extern int Nlocal_charge;
#define NDIM_MAX  3
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int Nwall;
extern int Pos_new_nodes;
extern int Plane_new_nodes;
#define NWALL_MAX_TYPE 50 
extern double Del_1[NWALL_MAX_TYPE];
extern double Size_x[NDIM_MAX];
extern double Energy;
extern int Lhard_surf;
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
int stencil_Njcomp_switch(int sten);
#define NSTEN        7
extern int Sten_Type[NSTEN];
extern int Nzone;
extern double **Vext_membrane;
extern int **Lsemiperm;
extern int **Xtest_reflect_TF;
extern int *Nwall_this_link;
extern int **Link_list;
extern int Nwall_type;
extern int *B2L_unknowns;
extern double *Deltac_b;
#define COULOMB      1
extern int Ipot_ff_c;
extern int L_Schur;
double gmin_double(double c);
double gmax_double(double c);
extern double *X_old;
extern int Nodes_old;
extern double *X2_old;
void post_process(double **x,char *output_file3,int *niters,double *time_save,int loop1,int binodal_flag);
int solve_problem(double **x,double **x2);
#define PICARD_NOX            3
#define NEWTON_NOX            1
#define NEWTON_BUILT_IN       0
int solve_problem_picard(double **x,double **x2);
#define PICARD_BUILT_IN       2
extern int NL_Solver;
extern int Lbinodal;
extern int Nnodes_box;
extern int Nunk_per_node;
extern int **Zero_density_TF;
void print_zeroTF(int **zero_TF,char *output_file);
extern double **Vext_static;
#define READ_VEXT_STATIC     3
extern int Restart_Vext;
extern double **Vext;
void print_vext(double **vext,char *output_file);
#define VERBOSE      3 
void setup_vext_coulomb_vol();
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define IN_WALL             -1
#define IN_BULK              0
#define REFLECT              2
#define PERIODIC             1
extern int Type_bc[NDIM_MAX][2];
extern int Vol_charge_flag;
extern double **Vext_coul;
void setup_integrals();
void setup_domain_multipliers();
void boundary_setup(char *output_file1);
void set_up_mesh(char *output_file1,char *output_file2);
extern int **Nel_hit2;
extern int **Nel_hit;
extern int *Comm_offset_unk;
extern int *Comm_offset_node;
extern int *Comm_unk_proc;
extern int *Comm_node_proc;
extern double *Area_IC;
void safe_free(void **ptr);
void safe_free(void **ptr);
extern int Lsteady_state;
void boundary_free(void);
void free_mesh_arrays(void);
void thermodynamics(char *output_file1);
void calc_stencils(void);
void calc_HS_diams();
#define NONE       -1
#define NONE      -1
#define NONE        -1
#define NONE        -1
extern int Type_func;
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Lmesh_refine;
extern double Esize_x[NDIM_MAX];
extern int Ndim;
void continuation_shift();
extern int Imain_loop;
#define NO_SCREEN    2 
extern int Iwrite;
extern int Nruns;
typedef struct Loca_Struct Loca_Struct;
extern struct Loca_Struct Loca;
void setup_polymer_cr();
#define N_NZCR_MAX   200   /* maximum # of non-zero's in direct correlation fn */
extern int Ncomp;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
extern double ***Rism_cr;
#define CMS_SCFT     1
#define CMS          0
extern int Type_poly;
void pot_parameters(char *output_file1);
extern int Mix_type;
void setup_nunk_per_node(char *output_file1);
void setup_stencil_logicals();
void read_input_file(char *input_file,char *output_file1);
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Num_Proc;
extern double Time_fill_av;
extern double Time_fill_first;
extern double Time_manager_av;
extern double Time_manager_first;
extern double Time_linsolver_av;
extern double Time_linsolver_first;
void dftmain(double *engptr);
