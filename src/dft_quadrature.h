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
#define PI    3.141592653589793238462643383279502884197169399375
void theta_midpoint(double **point,double *wt,int izone,int num_dim);
void get_radial_quadrature(double gauss_pt[],double gauss_wt[],int num_gp);
#define NSTEN        8
#define NZONE_MAX  10 
extern int Sten_Choice_R[NSTEN][NZONE_MAX];
#define THETA_CR_GENERAL_MSA  7
#define THETA_CR_RPM_MSA      3
#define THETA_PAIRPOT_RCUT    2
#define THETA_FN_R            1
void delta_midpoint(double **point,double *wt,int izone);
extern int MPsten_Npts_phi[NZONE_MAX];
extern int MPsten_Npts_arc[NZONE_MAX];
extern int MPsten_Npts_R[NZONE_MAX];
#define MIDPOINT_RULE 8
void delta_one_D_twelve(double **point,double *wt);
extern int Ndim;
#define ONE_D_TWELVE  7
void delta_seventy_two(double **point,double *wt);
#define SEVENTY_TWO   6
void delta_dodecahedron(double **point,double *wt);
#define DODECAHEDRON  5
void delta_icosahedron(double **point,double *wt);
#define ICOSAHEDRON   4
void delta_cube(double **point,double *wt);
#define CUBE          3
void delta_octahedron(double **point,double *wt);
#define OCTAHEDRON    2
void delta_tetrahedron(double **point,double *wt);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
#define TETRAHEDRON   1
extern int Sten_Choice_S[NSTEN][NZONE_MAX];
#define DELTA_FN_BOND         6
#define DELTA_FN_R            0
int get_integration_pts(int isten,int izone,double ***point_ptr,double **wt_ptr);
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
void set_gauss_quad(int ngp,double *gp,double *gw);
