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
extern double **Charge_x;
extern int Nlocal_charge;
#define NDIM_MAX  3
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int Pos_new_nodes;
int round_to_int(double x);
extern int Plane_new_nodes;
#define NWALL_MAX_TYPE 20 
extern double Del_1[NWALL_MAX_TYPE];
extern double Size_x[NDIM_MAX];
extern double Energy;
extern double *Vext_old;
extern double *X2_old;
extern double *X_old;
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
int stencil_Njcomp_switch(int sten);
#define NSTEN        8
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
extern double Time_InitGuess;
int gmax_int(int c);
int gmin_int(int c);
extern char *Outpath;
extern int L_Schur;
double gmin_double(double c);
double gmax_double(double c);
extern int Nodes_x[NDIM_MAX];
extern int Nodes_x_old[NDIM_MAX];
extern int Nnodes;
extern int Nodes_old;
#define FROM_MAIN 1
void post_process(double **x,int *niters,double *time_save,int loop1,int binodal_flag,int call_from_flag);
#define NEWTON_NOX            1
#define NEWTON_BUILT_IN       0
int solve_problem(double **x,double **x2);
void print_profile_box(double **x,char *outfile);
int solve_problem_picard(double **x,double **x2);
extern double NL_update_scalingParam;
#define PICNEWTON_BUILT_IN    4
#define PICNEWTON_NOX         5
#define PICARD_NOX            3
#define PICARD_BUILT_IN       2
extern int NL_Solver;
extern int Lbinodal;
extern int Nnodes_box;
extern int Nunk_per_node;
extern int **Zero_density_TF;
void print_zeroTF(int **zero_TF,char *ZeroTF_filename);
extern double **Vext_static;
#define READ_VEXT_STATIC     3
extern int Restart_Vext;
extern double **Vext;
void print_vext(double **vext,char *output_file);
#define FILES_DEBUG        2
#define FILES_EXTENDED     1 
#define SCREEN_VERBOSE     3 
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
extern int *List_coarse_nodes;
extern int Nnodes_coarse_loc;
extern int Nwall;
extern int *B2G_node_extra;
extern int Nnodes_box_extra;
void setup_global_surfaces();
extern int *L2G_node;
extern int Nnodes_per_proc;
void linsolver_setup_control();
extern void *LinProbMgr_manager;
void boundary_setup(char *file_echoinput);
void set_up_mesh(char *file_echoinput,char *output_file2);
#define NCOMP_MAX 5
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern int Grafted_TypeID[NCOMP_MAX];
extern int Grafted[NCOMP_MAX];
extern int Npol_comp;
extern int Grafted_Logical;
#define WJDC3        5 
extern int Lhard_surf;
extern int **Nel_hit2;
extern int **Nel_hit;
extern int *Comm_offset_unk;
extern int *Comm_offset_node;
extern int *Comm_unk_proc;
extern int *Comm_node_proc;
extern double *Area_IC;
void safe_free(void **ptr);
void safe_free(void **ptr);
#define DIFFUSIVE_INTERFACE 1
extern int Type_interface;
void boundary_free(void);
void free_mesh_arrays(void);
extern int Iwrite_files;
void thermodynamics(char *file_echoinput,int iwrite_screen,int iwrite_files);
void calc_stencils(void);
#define SCREEN_BASIC       1
#define SCREEN_ERRORS_ONLY  0 
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
void calc_HS_diams();
extern int Type_func;
typedef struct Loca_Struct Loca_Struct;
extern struct Loca_Struct Loca;
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Lmesh_refine;
extern double Esize_x[NDIM_MAX];
extern int Ndim;
void continuation_shift();
extern int Nruns;
extern int Imain_loop;
void setup_polymer_cr();
#define N_NZCR_MAX   200   /* maximum # of non-zero's in direct correlation fn */
extern int Ncomp;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern double ***Rism_cr;
#define CMS          0
extern int Type_poly;
void setup_pairPotentials(char *file_echoinput);
void setup_nunk_per_node(char *file_echoinput);
void setup_stencil_uattr_core_properties();
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_attr;
void setup_stencil_logicals();
#define FILENAME_LENGTH 300
extern char EchoInputFile_array[FILENAME_LENGTH];
void setup_params_for_dft(char *input_file,char *file_echoinput);
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Num_Proc;
extern double Time_NLSolve;
extern double Time_fill_av;
extern double Time_fill_first;
extern double Time_manager_av;
extern double Time_manager_first;
extern double Time_linsolver_av;
extern double Time_linsolver_first;
void dftmain(double *engptr);
