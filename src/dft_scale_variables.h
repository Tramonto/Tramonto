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
#define SCREEN_ERRORS_ONLY  0 
#define NCOMP_MAX 6
extern int Nmer[NCOMP_MAX];
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern double Rho_b[NCOMP_MAX];
#define NBLOCK_MAX   20 
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
extern int Npol_comp;
#define FILES_BASIC        0
#define SCREEN_NONE       -1 
void thermodynamics(char *file_echoinput,int iwrite_screen,int iwrite_files);
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern double Betamu[NCOMP_MAX];
#define NMER_MAX     200
extern double Betamu_chain[NMER_MAX];
#define WJDC3        5 
#define WJDC2        4 
#define WJDC         3
extern int Type_poly;
void calc_new_density(int icomp,char *file_echoinput);
#define NWALL_MAX_TYPE 20 
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
#define VEXT_HARD        1
extern int Ipot_wf_n[NWALL_MAX_TYPE];
extern int Nwall_type;
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_attr;
extern int Mix_type;
void scale_all_epsParams(double ratio);
void set_new_membrane_potential(double param_old,double param_new,int icomp);
void print_charge_surf(double **charge_w_sum,char *output_file);
void print_charge_vol(double *charge_els,char *output_file);
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define SCREEN_VERBOSE     3 
extern int Iwrite_screen;
#define NWALL_MAX 600 
extern double Elec_param_w[NWALL_MAX];
extern double *Charge_vol_els;
extern int Nelements_box;
extern int Vol_charge_flag;
extern double **Charge_w_sum_els;
extern int Surf_charge_flag;
void scale_elec_param(double ratio);
void scale_vext_epswf_terms(double ratio,int icomp,int iwall_type);
void sum_vext_epswf_terms();
void print_vext(double **vext,char *output_file);
#define VERBOSE      3 
extern int Iwrite;
extern double ***Vext_dash;
extern int Ndim;
extern int Nwall;
#define READ_VEXT_FALSE      0
extern int Lvext_dash;
extern double VEXT_MAX;
extern double **Vext_static;
extern double **Vext;
extern double ***Vext_perWallType;
extern int     WallType[NWALL_MAX]; /* array containing type number for each surface */
#define READ_VEXT_STATIC     3
extern int Restart_Vext;
extern int Nnodes_per_proc;
void scale_vext_temp(double ratio);
void calc_stencils(void);
extern int Lhard_surf;
void safe_free(void **ptr);
void safe_free(void **ptr);
typedef struct Stencil_Struct Stencil_Struct;
extern struct Stencil_Struct ***Stencil;
extern int Ncomp;
int stencil_Njcomp_switch(int sten);
#define NSTEN        8
extern int Sten_Type[NSTEN];
extern int Nzone;
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
void recalculate_stencils();
