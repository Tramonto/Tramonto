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
extern int Type_interface;
double calc_free_energy(FILE *fp,double **x);
void calc_force(FILE *fp,double **x,double fac_area);
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
void print_cont_variable(int cont_type,FILE *fp,int Loca_contID);
void print_cont_type(int cont_type,FILE *fp,int Loca_contID);
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
void print_gofr(char *output_file6,double *xold);
extern int Nlocal_charge;
extern int Nlink;
extern int Lprint_gofr;
void print_profile(char *output_file4,double *xold);
void collect_vext_old();
void collect_x_old(double **x,double *xold);
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
extern double *X2_old;
#define MINIMAL      0
#define NO_SCREEN    2 
extern int Iwrite;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define PRINT_RHO_0      0
extern int Print_rho_type;
#define FROM_LOCA 0
extern int Nruns;
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void post_process(double **x,int *niters,double *time_save,int loop1,int binodal_flag,int call_from_flag);
