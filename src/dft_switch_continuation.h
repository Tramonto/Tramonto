/* This file was automatically generated.  Do not edit! */
void print_cont_variable_user_plugin(int cont_type,FILE *fp,int Loca_contID);
void print_cont_variable_archived_plugin(int cont_type,FILE *fp,int Loca_contID);
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
#define NCOMP_MAX 5
extern double Charge_f[NCOMP_MAX];
#define PI    3.141592653589793238462643383279502884197169399375
#define NDIM_MAX  3
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
#define cyl3D_slit2D_pore               8
#define NWALL_MAX_TYPE 50 
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
#define atomic_centers                  3
#define finite_cyl_3D                   5
#define colloids_cyl_sphere             2
#define point_surface                   4
#define finite_planar_wall              1
#define smooth_planar_wall              0
typedef struct SurfaceGeom_Struct SurfaceGeom_Struct;
extern struct SurfaceGeom_Struct *SGeom;
extern int WallType[NWALL_MAX];
extern int Lwedge_cutout[NWALL_MAX];
extern int Lperiodic_overlay[NWALL_MAX_TYPE];
extern int Nperiodic_overlay[NWALL_MAX_TYPE];
#define NPERIODIC_MAX 4
extern int Llinear_overlay[NWALL_MAX_TYPE];
extern int Nlinear_overlay[NWALL_MAX_TYPE];
struct SurfaceGeom_Struct {
  int       surfaceTypeID;    /* ID of the type of surface */
  int       orientation;  /* orientation of the surface */
  double    *halfwidth;   /* planar surface params given in halfwidth */
  double    radius;       /* radius of spherical or cylindrical surface */
  double    halflength;   /* length of finite cylinders and pores */
  double    radius2;      /* a second radius for tapered pores or cylinders */
  int       Lwedge_cutout;    /* TRUE or FALSE for wedge cutout from basic surfac */
  double    angle_wedge_start;    /* angle as measured from x0 axis */
  double    angle_wedge_end;    /* angle as measured from x0 axis */
  int       Lrough_surface;    /* TRUE or FALSE for surface roughness */
  double    roughness;          /* maximum roughness amplitude */
  double    roughness_length;    /* lengthscale for the rougness */
  int       Lperiodic_overlay;    /* TRUE or FALSE for periodic function added to surface */
  int       Nperiodic_overlay;     /* The number of periodic functions to apply */
  double    orientation_periodic[NPERIODIC_MAX];    /* maximum amplitude for a cosine wave superimposed on a surface */
  double    amplitude[NPERIODIC_MAX];    /* maximum amplitude for a cosine wave superimposed on a surface */
  double    wavelength[NPERIODIC_MAX];    /* desired wavelength of cosine wave superimposed on a surface */
  double    origin_PeriodicFunc[NPERIODIC_MAX];     /* The origin of periodic functions to apply */
  int       Llinear_overlay;    /* TRUE or FALSE for linear function added to surface */
  int       Nlinear_overlay;     /* The number of linear functions to apply */
  double    orientation_linear[NPERIODIC_MAX];    /* maximum amplitude for a linear function superimposed on a surface */
  double    slope[NPERIODIC_MAX];    /* maximum amplitude for a linear function superimposed on a surface */
  double    origin_LinearFunc[NPERIODIC_MAX];     /* The origin of linear functions to apply */
  double    endpoint_LinearFunc[NPERIODIC_MAX];     /* The end point of linear functions to apply */
  int    *ReflectionsAreIndependent;  /* TRUE or FALSE for treating special boundary conditions */
};
double print_cont_variable(int cont_type,FILE *fp,int Loca_contID);
void print_cont_type_user_plugin(int cont_type,FILE *fp,int Loca_contID);
void print_cont_type_archived_plugin(int cont_type,FILE *fp,int Loca_contID);
#define SWITCH_MU    3
#define SWITCH_ION   2
#define SWITCH_ALLTYPES_ICOMP 1
#define SWITCH_RHO   0
extern int Npol_comp;
#define SWITCH_BULK_OUTPUT 5
#define SWITCH_ALLTYPES 4
extern int Print_rho_switch;
extern int Type_attr;
extern int Ndim;
#define REFLECT              2
extern int Type_bc[NDIM_MAX][2];
#define SWITCH_SURFACE_SEP   0
extern int Print_mesh_switch;
void print_cont_type(int cont_type,FILE *fp,int Loca_contID);
void thermodynamics(char *output_file1);
void recalculate_stencils();
void setup_polymer_cr();
void calc_InvR_params();
void calc_HS_diams();
void set_new_membrane_potential(double param_old,double param_new,int icomp);
void setup_wall_wall_potentials();
extern int Lprint_pmf;
void scale_vext_temp(double ratio);
void scale_vext_epswf(double ratio,int icomp,int iwall);
void scale_elec_param(double ratio);
void setup_pairPotentials(char *output_file1);
void scale_all_epsParams(double ratio);
extern int Ncomp;
#define BH_DIAM             1
extern int Type_hsdiam;
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
extern int Nwall;
#define CMS          0
extern int Type_func;
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern double Temp_elec;
#define COULOMB      1
extern int Ipot_ff_c;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
void adjust_dep_params(int cont_type,int Loca_contID,double param_old,double param_new,char *output_file1);
void assign_param_user_plugin(int cont_type,int Loca_contID,double param);
void assign_param_archived_plugin(int cont_type,int Loca_contID,double param);
#define NMER_MAX     200
extern int Unk2Comp[NMER_MAX];
extern int Ntype_mer;
void assign_parameter_tramonto(int cont_type,double param,int Loca_contID);
double get_init_param_user_plugin(int cont_type,int Loca_contID);
double get_init_param_archived_plugin(int cont_type,int Loca_contID);
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
#define CONT_SIGMAFF_IJ    10   /* Fluid-Fluid Interaction Diameter for IJ term */
extern double **Vext_membrane;
#define CONT_SEMIPERM_IJ   9  /* Vext_membrane */
#define CONT_ELECPARAM_ALL 8  /* Charged surface params */
extern double Elec_param_w[NWALL_MAX];
#define CONT_ELECPARAM_I   7  /* Charged surface params */
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
#define CONT_EPSFF_IJ      6   /* Fluid-Fluid Energy Params for IJ term */
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
#define CONT_EPSWF_IJ      5    /* Wall-Fluid Energy Params for IJ term */
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Eps_w[NWALL_MAX_TYPE];
extern int Mix_type;
#define CONT_EPSW_I        4    /* Wall-Wall Energy Params for wall I */
extern double Betamu[NCOMP_MAX];
extern double Betamu_chain[NMER_MAX];
#define WJDC3        5 
#define WJDC2        4 
#define WJDC         3
#define CONT_BETAMU_I      3  /* Vary chemical potential for species I */
extern double Rho_seg_b[NMER_MAX];
extern int SegAll_to_Poly[NMER_MAX];
extern int Nseg_tot;
#define NCONT_MAX          2 /* the maximum number of solutions possible for use with Loca */
extern int Cont_ID[NCONT_MAX][2];
extern double Rho_b[NCOMP_MAX];
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_poly;
#define CONT_RHO_I         2
extern double Temp;
#define CONT_TEMP          1   /* State Parameters */
extern int Plane_new_nodes;
extern double Size_x[NDIM_MAX];
#define CONT_MESH          0   /* mesh size */
#if defined(LOCA_MF) && !defined(_CON_CONST_H_)
#include <mf.h>
#endif
double get_init_param_value(int cont_type,int);
#if !defined(_CON_CONST_H_)
double get_init_param_value(int cont_type,int);
#endif
double get_init_param_value(int cont_type,int Loca_contID);
