/* This file was automatically generated.  Do not edit! */
void thermodynamics(char *output_file1);
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
typedef struct Loca_Struct Loca_Struct;
struct Loca_Struct {
  int    method;      /* Continuation method                          */
  int    cont_type1;  /* flag specifying the continuation parameter   */
  int    cont_type2;  /* flag specifying the second (free) parameter  */
  int    num_steps;   /* maximum number of continuation steps to take */
  double aggr;        /* step size control parameter                  */
  double step_size;   /* initial continuation step size               */
};
extern struct Loca_Struct Loca;
extern int **Zero_density_TF;
void print_zeroTF(int **zero_TF,char *output_file);
void print_vext(double **vext,char *output_file);
#define VERBOSE      3 
extern int Iwrite;
void boundary_setup(char *output_file1);
void set_up_mesh(char *output_file1,char *output_file2);
extern int *Comm_offset_unk;
extern int *Comm_offset_node;
extern int *Comm_unk_proc;
extern int *Comm_node_proc;
extern double *Area_IC;
void safe_free(void **ptr);
void safe_free(void **ptr);
extern int Ndim;
#define DIFFUSIVE_INTERFACE 1
extern int Type_interface;
void boundary_free(void);
void free_mesh_arrays(void);
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern double **Vext;
extern int Nnodes_per_proc;
void scale_elec_param(double ratio);
void scale_vext_epswf(double ratio,int icomp,int iwall);
#define NMER_MAX     100
extern int Unk2Comp[NMER_MAX];
extern int Nseg_tot;
void scale_vext_temp(double ratio);
extern int Nwall;
void recalculate_stencils();
void setup_polymer_cr();
#define CMS          0
void calc_InvR_params();
void calc_HS_diams();
extern int Type_func;
void pot_parameters(char *output_file1);
#define VEXT_HARD        1
#define NWALL_MAX_TYPE 50 
extern int Ipot_wf_n[NWALL_MAX_TYPE];
#define LJ12_6       2
extern int Ipot_ff_n;
extern double Temp_elec;
#define COULOMB      1
extern int Ipot_ff_c;
#define NCOMP_MAX 5
void assign_parameter_tramonto(int cont_type,double param);
extern double Crfac;
#define CONT_CRFAC  19  /* continuous mixing of two cr files */
extern int WallType[NWALL_MAX];
extern double WallParam[NWALL_MAX_TYPE];
#define CONT_WALLPARAM  18  /* Vext_membrane */
extern double **Vext_membrane;
#define CONT_SEMIPERM   17  /* Vext_membrane */
#define CONT_SCALE_CHG   16  /* Charged surface params */
#define CONT_SCALE_EPSFF 15
#define CONT_EPSFF_ALL   14   
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
#define CONT_EPSFF_00    13   /* Fluid-Fluid Energy Params */
#define CONT_SCALE_EPSWF 12
#define CONT_EPSWF_ALL_0 11 
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
#define CONT_EPSWF00     10    /* Wall-Fluid Energy Params */
#define CONT_SCALE_EPSW  9
extern int Nwall_type;
#define CONT_EPSW_ALL    8
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Eps_w[NWALL_MAX_TYPE];
extern int Mix_type;
#define CONT_EPSW_0      7    /* Wall-Wall Energy Params */
#define CONT_BETAMU_1 21  /* Vary chemical potential for species 1 */
extern double Betamu[NCOMP_MAX];
extern double Betamu_chain[NMER_MAX];
#define WJDC3        5 
#define WJDC2        4 
#define WJDC         3
#define CONT_BETAMU_0 20  /* Vary chemical potential for species 0 */
extern double Scale_fac;
#define CONT_SCALE_RHO   6
#define CONT_LOG_RHO_ALL 5 
#define CONT_LOG_RHO_0   4 
extern double Rho_seg_b[NMER_MAX];
extern int Npol_comp;
extern int Ncomp;
#define NONE       -1
#define NONE      -1
#define NONE        -1
#define NONE        -1
extern int Type_poly;
#define CONT_RHO_ALL     3
extern double Rho_b[NCOMP_MAX];
#define CONT_RHO_0       2
extern double Temp;
#define CONT_TEMP        1   /* State Parameters */
#define CONT_MESH        0   /* mesh size */
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if defined(LOCA_MF) && !defined(_CON_CONST_H_)
#include <mf.h>
#endif
double get_init_param_value(int cont_type);
#if !defined(_CON_CONST_H_)
double get_init_param_value(int cont_type);
#endif
double get_init_param_value(int cont_type);
