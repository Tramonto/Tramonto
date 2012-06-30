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
#define NZONE_MAX  10 
extern int Num_Proc;
#define NWALL_MAX_TYPE 20 
#define NWALL_MAX 600 
void error_check(void);
#define REFLECT              2
#define LAST_NODE            3
#define IN_BULK              0
#define CONT_BETAMU_I_NEW  11 /* Vary chemical potential for species I...holding densities of other species constant */
#define CONT_BETAMU_I      3  /* Vary chemical potential for species I */
#define CONT_EPSWF_ALL	        105
#define CONT_EPSWF_IJ      5    /* Wall-Fluid Energy Params for IJ term */
#define CONT_LOG_RHO_I          100
#define CONT_RHO_I         2
extern int NID_Cont;
#define NCONT_MAX          2 /* the maximum number of solutions possible for use with Loca */
extern int Cont_ID[NCONT_MAX][2];
extern int Lbinodal;
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
extern int Pos_new_nodes;
extern int Plane_new_nodes;
extern double Del_1[NWALL_MAX_TYPE];
extern int Nruns;
extern double Az_tolerance;
extern int Max_gmres_iter;
extern double Az_ilut_fill_param;
extern int Az_preconditioner;
extern int Az_scaling;
extern int Az_kspace;
extern int Az_solver;
extern int L_Schur;
extern int Load_Bal_Flag;
extern double NL_abs_tol_picard,NL_rel_tol_picard;
extern double NL_update_scalingParam;
extern double NL_abs_tol,NL_rel_tol;
#define CALC_ALL_FIELDS   1
#define PICARD_BUILT_IN       2
#define NCOMP_MAX 5
extern double Scale_fac_WJDC[NCOMP_MAX][NCOMP_MAX];
#define WJDC2        4 
#define WJDC         3
extern int Analyt_WJDC_Jac;
extern int ATTInA22Block;
extern int Physics_scaling;
extern int Max_NL_iter;
extern int NL_Solver;
extern int Dim_1Dbc;
extern double X_1D_bc;
extern double Jac_threshold;
extern int Lcut_jac;
extern double Jac_grid;
extern int Coarser_jac;
#define SCREEN_ERRORS_ONLY  0 
extern int Mesh_coarsening;
extern double Rmax_zone[5];
extern int Nzone;
#define FILES_DEBUG_MATRIX 3 
#define VERBOSE_MATRIX    5 
#define SCREEN_NONE       -1 
#define NO_SCREEN    4 
#define FILES_DEBUG        2
#define SCREEN_VERBOSE     3 
#define VERBOSE      3 
#define FILES_EXTENDED     1 
#define EXTENDED     2
#define FILES_BASIC        0
extern int Iwrite_files;
#define SCREEN_BASIC       1
extern int Iwrite_screen;
#define DENSITIES    1 
#define MINIMAL      0
extern int Iwrite;
extern int Print_mesh_switch;
extern int Print_rho_switch;
extern int Print_rho_type;
extern int Lprint_pmf;
extern int Lprint_gofr;
extern int Lcount_reflect;
extern int Lper_area;
extern double Rho_max;
extern char *Vext_filename2;
#define FILENAME_LENGTH 300
extern char vext_file2_array[FILENAME_LENGTH];
#define READ_VEXT_STATIC     3
#define READ_VEXT_SUMTWO     2
extern char *Vext_filename;
extern char vext_file_array[FILENAME_LENGTH];
#define READ_VEXT_FALSE      0
extern int Restart_Vext;
extern int Nmissing_densities;
#define RESTART_FEWERCOMP  4
extern int Restart;
#define NSTEPS_MAX 10
extern double Rho_step[NCOMP_MAX][NSTEPS_MAX];
extern double Xend_step[NSTEPS_MAX];
extern double Xstart_step[NSTEPS_MAX];
extern int Orientation_step[NSTEPS_MAX];
#define CHOP_RHO_STEP    4
#define CHOP_RHO         3
extern int Nsteps;
#define STEP_PROFILE     2
extern int Iguess_fields;
extern int Iguess;
extern double *Lseg_IC;
extern double *Pore_rad_R_IC;
extern double *Pore_rad_L_IC;
extern int Nseg_IC;
extern int Geom_flag;
extern double Velocity;
extern double D_coef[NCOMP_MAX];
#define DIFFUSIVE_INTERFACE 1
extern int Linear_transport;
extern int L1D_bc;
extern double *Dielec_wall;
extern double Dielec_X;
extern double Dielec_pore;
extern double Dielec_bulk;
#define DIELEC_WF_PORE     2
extern double DielecConst_plasma;
extern double Temp_K_plasma;
extern double Sigma_Angstroms_plasma;
extern int Type_dielec;
#define POINT_CHARGE 0
extern int Charge_type_local;
extern int Charge_type_atoms;
#define ATOMIC_CHARGE    3
extern double **Charge_x;
extern double *Charge_Diam;
extern double *Charge;
extern int Nlocal_charge;
extern int Type_bc_elec[NWALL_MAX_TYPE];
extern int Ipot_wf_c;
extern double X_const_mu;
extern double Elec_pot_RTF;
extern double Elec_pot_LBB;
#define COULOMB      1
extern double Rho_b_RTF[NCOMP_MAX];
extern double Rho_b_LBB[NCOMP_MAX];
extern double Rho_b[NCOMP_MAX];
extern int Lconstrain_interface;
extern int Grad_dim;
extern double **Vext_membrane;
extern int **Lsemiperm;
extern double Cr_rad_hs[NCOMP_MAX][NCOMP_MAX];
extern double Cr_break[2];
extern char *Cr_file2;
extern char cr_file2_array[FILENAME_LENGTH];
extern char *Cr_file;
extern char cr_file_array[FILENAME_LENGTH];
extern double Crfac;
extern int Ncr_files;
#define CMS          0
void setup_chain_architecture(char *poly_file,FILE *fpecho);
extern char *Poly_file_name;
extern char poly_file_array[FILENAME_LENGTH];
#define SCFT         6	
#define POLY_ARCH_FILE 0
extern int Type_poly_arch;
extern double Rho_g[NCOMP_MAX];
extern int Graft_wall[NCOMP_MAX];
extern int GraftedWall_TF[NWALL_MAX_TYPE];
extern int Grafted_TypeID[NCOMP_MAX];
extern int Grafted_SegIDAll[NCOMP_MAX];
extern int Grafted_SegID[NCOMP_MAX];
extern int Grafted[NCOMP_MAX];
#define NBLOCK_MAX   20 
extern int Poly_to_Type[NCOMP_MAX][NBLOCK_MAX];
#define NMER_MAX     200
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern int SegType_per_block[NCOMP_MAX][NBLOCK_MAX];
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
extern int Poly_to_Ntype[NCOMP_MAX];
extern int Type_mer_to_Pol[NBLOCK_MAX];
extern int Nmer_t_total[NBLOCK_MAX];
extern int Ntype_mer;
extern int Nseg_per_block[NCOMP_MAX][NBLOCK_MAX];
extern int Nmer[NCOMP_MAX];
extern int Nblock[NCOMP_MAX];
extern int Npol_comp;
extern int Grafted_Logical;
extern int Unk2Comp[NMER_MAX];
extern double YukawaK_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double EpsYukawa_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Cut_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double YukawaK_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double EpsYukawa_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Cut_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
#define atomic_centers                  3
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Rho_w[NWALL_MAX_TYPE];
extern double Npow_ff[NCOMP_MAX][NCOMP_MAX];
extern double YukawaK_ff[NCOMP_MAX][NCOMP_MAX];
extern double EpsYukawa_ff[NCOMP_MAX][NCOMP_MAX];
#define PAIR_rNandYUKAWA_CS   9
#define PAIR_r18andYUKAWA_CS  8
#define PAIR_r12andYUKAWA_CS  7
#define PAIR_LJandYUKAWA_CS   6
#define PAIR_EXP_CS	      4
#define PAIR_YUKAWA_CS        3
extern double Bond_ff[NCOMP_MAX][NCOMP_MAX];
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern int Lpolarize[NCOMP_MAX];
extern double Pol[NCOMP_MAX];
extern double Charge_f[NCOMP_MAX];
extern double Mass[NCOMP_MAX];
extern double HS_diam[NCOMP_MAX];
#define MANUAL_HS_DIAM         2
#define VEXT_HARD        1
#define VEXT_NONE          0
extern int Mix_type;
extern int Ncomp;
extern int Type_uwwPot;
extern int Ipot_ww_n[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern int Vext_PotentialID[NWALL_MAX_TYPE];
extern int Type_vext[NWALL_MAX_TYPE];
extern int Lhard_surf;
extern int Ipot_wf_n[NWALL_MAX_TYPE];
extern double Angle_wedge_end[NWALL_MAX];
extern double Angle_wedge_start[NWALL_MAX];
extern int Lwedge_cutout[NWALL_MAX];
extern int read_wedge;
#define NPERIODIC_MAX 4
extern double EndpointLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern double OriginLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern double SlopeLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern int OrientationLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern int Nlinear_overlay[NWALL_MAX_TYPE];
extern int Llinear_overlay[NWALL_MAX_TYPE];
extern int read_linear;
extern double OriginPeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern double WavelengthPeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern double AmplitudePeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern int OrientationPeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern int Nperiodic_overlay[NWALL_MAX_TYPE];
extern int Lperiodic_overlay[NWALL_MAX_TYPE];
extern int read_periodic;
extern double Rough_length[NWALL_MAX_TYPE];
extern double Rough_param_max[NWALL_MAX_TYPE];
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Lrough_surf[NWALL_MAX_TYPE];
extern int read_rough;
extern int Lapply_offset[3];
extern double WallParam_3[NWALL_MAX_TYPE];
extern double WallParam_2[NWALL_MAX_TYPE];
#define point_surface                   4
extern double WallParam[NWALL_MAX_TYPE];
extern int Orientation[NWALL_MAX_TYPE];
extern int Surface_type[NWALL_MAX_TYPE];
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern int **Xtest_reflect_TF;
extern int Lauto_size;
extern int Lauto_center;
extern int Nlink;
extern int Nwall;
extern int Nwall_type;
extern int Ipot_ff_c;
#define PAIR_COULOMB          2
#define PAIR_COULOMB_CS       1
#define LJ12_6       2
#define HARD_SPHERE  1
#define IDEAL_GAS    0
extern int Ipot_ff_n;
#define WTC          2
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
#define WJDC3        5 
extern int Type_poly;
extern int Type_coul;
extern int Type_pairPot;
extern int Type_attr;
extern int Type_hsdiam;
extern int Type_func;
#define PERIODIC             1
#define NDIM_MAX  3
extern int Type_bc[NDIM_MAX][2];
extern double Esize_x[NDIM_MAX];
extern double Size_x[NDIM_MAX];
extern int Ndim;
extern double VEXT_MAX;
extern double Dielec_ref;
extern double Temp;
extern double Density_ref;
extern double Length_ref;
void read_junk(FILE *fpinput,FILE *fpecho);
#define UNIFORM_INTERFACE  0
extern int Type_interface;
extern int LBulk;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int LDeBroglie;
void read_input_file(FILE *fpinput,FILE *fpecho);
