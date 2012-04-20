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
#define cyl3D_slit2D_pore               8
#define cyl2D_sphere3D_pore             7
#define atomic_centers                  3
#define finite_cyl_3D                   5
#define colloids_cyl_sphere             2
#define point_surface                   4
#define finite_planar_wall              1
#define NWALL_MAX_TYPE 20 
#define NPERIODIC_MAX 4
extern double EndpointLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern double OriginLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern double SlopeLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern int OrientationLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern int Nlinear_overlay[NWALL_MAX_TYPE];
extern double OriginPeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern double WavelengthPeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern double AmplitudePeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern int OrientationPeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern int Nperiodic_overlay[NWALL_MAX_TYPE];
#define NWALL_MAX 600 
extern double Angle_wedge_end[NWALL_MAX];
extern double Angle_wedge_start[NWALL_MAX];
extern double Rough_length[NWALL_MAX_TYPE];
extern int Lrough_surf[NWALL_MAX_TYPE];
extern int Orientation[NWALL_MAX_TYPE];
#define smooth_planar_wall              0
extern int Lwedge_cutout[NWALL_MAX];
extern int Llinear_overlay[NWALL_MAX_TYPE];
extern int Lperiodic_overlay[NWALL_MAX_TYPE];
extern int Surface_type[NWALL_MAX_TYPE];
extern double *Poly_graft_dist;
#define CONT_LOG_RHO_I          100
#define CONT_RHO_I         2
#define NCOMP_MAX 5
extern double Rho_b_RTF[NCOMP_MAX];
extern double Rho_b_LBB[NCOMP_MAX];
extern int Ntype_mer;
extern double Rho_w[NWALL_MAX_TYPE];
extern double *Dielec_wall;
extern double Dielec_pore;
extern double Dielec_bulk;
extern double Elec_pot_RTF;
extern double Elec_pot_LBB;
#define CONT_EPSFF_ALL		107
#define CONT_EPSFF_IJ      6   /* Fluid-Fluid Energy Params for IJ term */
#define CONT_EPSWF_ALL	        105
#define CONT_EPSWF_IJ      5    /* Wall-Fluid Energy Params for IJ term */
#define CONT_EPSW_ALL		104
#define CONT_EPSW_I        4    /* Wall-Wall Energy Params for wall I */
#define CONT_TEMP          1   /* State Parameters */
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
#define USE_LOCA
extern double **Vext_membrane;
extern double EpsYukawa_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double EpsYukawa_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double EpsYukawa_ff[NCOMP_MAX][NCOMP_MAX];
extern double X_1D_bc;
extern double Rmax_zone[5];
extern double Dielec_X;
#define DIELEC_WF_PORE     2
extern int Type_dielec;
extern double Elec_param_w[NWALL_MAX];
extern int WallType[NWALL_MAX];
extern int Type_bc_elec[NWALL_MAX_TYPE];
extern double X_const_mu;
#define UNIFORM_INTERFACE  0
extern int Type_interface;
extern double YukawaK_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Cut_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern int Mix_type;
extern double YukawaK_ff[NCOMP_MAX][NCOMP_MAX];
#define PAIR_rNandYUKAWA_CS   9
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
extern double Bond_ff[NCOMP_MAX][NCOMP_MAX];
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern double YukawaK_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
#define PAIR_r18andYUKAWA_CS  8
#define PAIR_r12andYUKAWA_CS  7
#define PAIR_LJandYUKAWA_CS   6
#define PAIR_EXP_CS	      4
#define PAIR_YUKAWA_CS        3
extern int Type_uwwPot;
extern double Cut_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double WallParam_3[NWALL_MAX_TYPE];
extern double WallParam_2[NWALL_MAX_TYPE];
extern double WallParam[NWALL_MAX_TYPE];
#define NDIM_MAX  3
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern double Esize_x[NDIM_MAX];
extern double Size_x[NDIM_MAX];
extern char *DensityFile2;
#define FILENAME_LENGTH 300
extern char DensityFile2_array[FILENAME_LENGTH];
extern int Lbinodal;
extern char *DensityFile;
extern char DensityFile_array[FILENAME_LENGTH];
extern int Lprint_scaleFacWJDC;
extern int Physics_scaling;
extern int Nzone;
#define JAC_ZONES_SETFIXED_ESIZE       5
extern int Coarser_jac;
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
#define HARD_SPHERE  1
extern int Ipot_ff_n;
extern int Nnodes_per_el_S;
extern int Ndim;
extern int Nnodes_per_el_V;
#define SCREEN_ERRORS_ONLY  0 
#define SCREEN_NONE       -1 
extern double Sigma_Angstroms_plasma;
#define EPSILON_0  8.85419e-12  /* C^2 J^-1 m^-1 */
extern double DielecConst_plasma;
extern double Temp_K_plasma;
#define PI    3.141592653589793238462643383279502884197169399375
extern double Temp_elec;
extern double Rho_b[NCOMP_MAX];
extern int Ncomp;
extern double Rho_t;
#define CMS_SCFT     1
extern double Rough_param_max[NWALL_MAX_TYPE];
#define MAX_ROUGH_BLOCK 100
extern double Rough_precalc[NWALL_MAX_TYPE][MAX_ROUGH_BLOCK][MAX_ROUGH_BLOCK];
#if defined(DEC_ALPHA)
#define POW_INT powii
#endif
#if !(defined(DEC_ALPHA))
#define POW_INT (int)pow
#endif
extern int read_rough;
extern int Link[NWALL_MAX];
extern int *Nwall_this_link;
extern int Nwall;
extern int Nlink;
extern int **Link_list;
extern int Lmesh_refine;
#define SCREEN_VERBOSE     3 
extern int Iwrite_screen;
void setup_other_run_constants();
void fill_surfGeom_struct();
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
typedef struct SurfaceGeom_Struct SurfaceGeom_Struct;
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
extern struct SurfaceGeom_Struct *SGeom;
extern int Nwall_type;
extern int ***pol_sym_tmp;
void safe_free(void **ptr);
void safe_free(void **ptr);
void setup_chain_indexing_arrays(int nseg,int nmer_max,FILE *fpecho);
extern int Nmer[NCOMP_MAX];
extern int Npol_comp;
void make_density_params_dimensionless();
extern double Density_ref;
void make_dielecConst_params_dimensionless();
extern double Dielec_ref;
void make_elecPot_params_dimensionless();
#define E_CONST    1.60219e-19  /* C */
#define KBOLTZ     1.3807e-23  /* J K-1 */
extern double Potential_ref;
#define PAIR_COULOMB          2
#define PAIR_COULOMB_CS       1
extern int Type_pairPot;
void make_energy_params_dimensionless();
void make_length_params_dimensionless();
extern double Length_ref;
void broadcast_input();
extern int Num_Proc;
extern int Flag_mV_elecpot;
extern int Type_coul;
extern double Temp;
extern char *Poly_file_name;
void setup_chain_architecture(char *poly_file,FILE *fpecho);
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_poly;
void read_input_file(FILE *fpinput,FILE *fpecho);
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Open_GUI;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void setup_params_for_dft(char *input_file,char *file_echoinput);
