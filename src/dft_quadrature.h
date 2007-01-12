/* This file was automatically generated.  Do not edit! */
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
#define PI    M_PI
void theta_midpoint(double **point,double *wt,int izone,int num_dim);
void get_radial_quadrature(double gauss_pt[],double gauss_wt[],int num_gp);
#define NSTEN        8
#define NZONE_MAX  10 
extern int Sten_Choice_R[NSTEN][NZONE_MAX];
#define THETA_CHARGE   3
#define U_ATTRACT      2
#define THETA_FN       1
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
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
#define TETRAHEDRON   1
extern int Sten_Choice_S[NSTEN][NZONE_MAX];
#define DELTA_FN       0
int get_integration_pts(int isten,int izone,double ***point_ptr,double **wt_ptr);
void set_gauss_quad(int ngp,double *gp,double *gw);
