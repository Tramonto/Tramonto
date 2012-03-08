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
extern double *X_old;
void calc_flux(FILE *fp,char *output_flux,double *X_old);
#define NCOMP_MAX 5
extern double Betamu[NCOMP_MAX];
double calc_free_energy(FILE *fp,double **x);
void calc_force(FILE *fp,double **x,double fac_area);
#define DIFFUSIVE_INTERFACE 1
extern int Type_interface;
void calc_fluid_charge(FILE *fp,double **x);
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_coul;
void calc_adsorption(FILE *fp,double **x);
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int Lcount_reflect;
#define REFLECT              2
extern int Type_bc[NDIM_MAX][2];
extern int Ndim;
void setup_domain_multipliers();
extern double Scale_fac_WJDC[NCOMP_MAX][NCOMP_MAX];
extern int **Nseg_type_pol;
extern int Npol_comp;
extern int Lprint_scaleFacWJDC;
#define CONT_MESH          0   /* mesh size */
double print_cont_variable(int cont_type,FILE *fp,int Loca_contID);
void print_cont_type(int cont_type,FILE *fp,int Loca_contID);
#define NCONT_MAX          2 /* the maximum number of solutions possible for use with Loca */
extern int Cont_ID[NCONT_MAX][2];
#define CONT_BETAMU_I_NEW  11 /* Vary chemical potential for species I...holding densities of other species constant */
#define CONT_RHO_I         2
#define SWITCH_NO_STATEOUT -1
extern int Print_rho_switch;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define CONT_BETAMU_I      3  /* Vary chemical potential for species I */
extern int Nwall;
#define FROM_MAIN 1
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
void print_gofr(char *GofR_Filename,double *xold);
extern int Nlocal_charge;
extern int Nlink;
extern int Lprint_gofr;
void print_profile(char *Density_FileName,double *xold);
void collect_vext_old();
void collect_x_old(double **x,double *xold);
extern int Ncomp;
extern double *Vext_old;
extern int Nunk_per_node;
extern int Nnodes;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern double *X2_old;
#define SCREEN_ERRORS_ONLY  0 
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
#define PRINT_RHO_0      0
extern int Print_rho_type;
extern int Imain_loop;
extern int Nruns;
extern char *OutputFileDir;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void post_process(double **x,int *niters,double *time_save,int loop1,int binodal_flag,int call_from_flag);
