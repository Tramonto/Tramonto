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
#define CONT_WALLPARAM  18  /* Vext_membrane */
extern double **Vext_membrane;
#define CONT_SEMIPERM   17  /* Vext_membrane */
#define CONT_SCALE_CHG   16  /* Charged surface params */
#define CONT_EPSFF_ALL   14   
#define CONT_EPSFF_00    13   /* Fluid-Fluid Energy Params */
#define CONT_SCALE_EPSFF 15
#define CONT_EPSWF_ALL_0 11 
#define CONT_EPSWF00     10    /* Wall-Fluid Energy Params */
#define NCOMP_MAX 5
#define NWALL_MAX_TYPE 50 
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
#define CONT_SCALE_EPSWF 12
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Eps_w[NWALL_MAX_TYPE];
extern int Mix_type;
#define CONT_EPSW_ALL    8
#define CONT_EPSW_0      7    /* Wall-Wall Energy Params */
#define CONT_SCALE_EPSW  9
#define CONT_BETAMU_1 21  /* continuous mixing of two cr files */
#define NMER_MAX     100
extern double Betamu_chain[NMER_MAX];
#define WJDC         3
extern int Type_poly;
#define CONT_BETAMU_0 20  /* continuous mixing of two cr files */
extern double Betamu[NCOMP_MAX];
#define SWITCH_MU    3
extern double Charge_f[NCOMP_MAX];
#define PI    M_PI
#define COULOMB      1
#define SWITCH_ION   2
extern double P_over_po;
#define SWITCH_RELP  1
extern int Print_rho_switch;
extern double Rho_b[NCOMP_MAX];
#define CONT_LOG_RHO_ALL 5 
#define CONT_LOG_RHO_0   4 
#define CONT_RHO_ALL     3
#define CONT_RHO_0       2
extern double Scale_fac;
#define CONT_SCALE_RHO   6
extern double Crfac;
#define CONT_CRFAC  19  /* continuous mixing of two cr files */
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
extern int Type_attr;
extern double Temp_elec;
extern double Temp;
extern int Ipot_ff_c;
#define CONT_TEMP        1   /* State Parameters */
extern double WallParam[NWALL_MAX_TYPE];
extern int Plane_new_nodes;
#define NWALL_MAX 600 
extern int WallType[NWALL_MAX];
extern int Nwall;
#define SWITCH_SURFACE_SEP   0
extern int Print_mesh_switch;
#define CONT_MESH        0   /* mesh size */
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern double *X_old;
void calc_flux(FILE *fp,char *output_flux,double *X_old);
#define DIFFUSIVE_INTERFACE 1
extern int Lsteady_state;
double calc_free_energy(FILE *fp,double **x);
void calc_force(FILE *fp,double **x,double fac_area);
void calc_fluid_charge(FILE *fp,double **x);
#define NONE       -1
#define NONE      -1
#define NONE        -1
#define NONE        -1
extern int Type_coul;
void calc_adsorption(FILE *fp,double **x);
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int Lcount_reflect;
#define REFLECT              2
extern int Type_bc[NDIM_MAX][2];
extern int Ndim;
void setup_domain_multipliers();
extern int Nruns;
void print_cont_variable(int cont_type,FILE *fp);
void print_cont_type(int cont_type,FILE *fp);
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
void safe_free(void **ptr);
void safe_free(void **ptr);
void print_gofr(char *output_file6);
extern int Nlink;
extern int Lprint_gofr;
void print_profile(char *output_file4);
#define MINIMAL      0
void collect_vext_old();
void collect_x_old(double **x);
extern int Ncomp;
extern double *Vext_old;
extern int Nunk_per_node;
extern int Nnodes;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
#define NO_SCREEN    2 
extern int Iwrite;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define PRINT_RHO_0      0
extern int Print_rho_type;
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void post_process(double **x,char *output_file3,int *niters,double *time_save,int loop1,int binodal_flag);
