/* This file was automatically generated.  Do not edit! */
void print_charge_surf(double **charge_w_sum,char *output_file);
void print_charge_vol(double *charge_els,char *output_file);
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
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define NWALL_MAX 600 
extern double Elec_param_w[NWALL_MAX];
extern double *Charge_vol_els;
extern int Nelements_box;
extern int Vol_charge_flag;
extern double **Charge_w_sum_els;
extern int Surf_charge_flag;
void scale_elec_param(double ratio);
void scale_vext_epswf(double ratio,int icomp,int iwall);
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
