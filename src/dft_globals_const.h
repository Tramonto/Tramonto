/*
//@HEADER
// ********************************************************************
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/*
 *  File:  dft_globals_const.h
 *
 *  This file contains the extern statements of some of the key global 
 *  varibles that are used often thoughout the code. The file 
 *  "dft_globals.h" should be included in the main program file, while
 *  this file should be included in any other files that use one of 
 *  these global variables.
 *
 *  It also includes some standard "C" libraries that need to be
 *  included for common things such as print statements, square-roots.
 */
/****************************************************************************/

/*** Standard "C" include files ***/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef HAS_VALUES_H
#include <values.h> /* for PI called (M_PI) */
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

/*#ifdef __cplusplus
extern "C" {
#endif*/

/****************************************************************************/
/* Machine specific definitions */

/* Integer exponentialtion requires powii on a DEC */
#ifdef DEC_ALPHA
#define POW_INT powii
#else
#define POW_INT (int)pow
#endif

/* Integer exponentialtion requires powii on a DEC */
#ifdef DEC_ALPHA
#define POW_DOUBLE_INT powi
#else
#define POW_DOUBLE_INT pow
#endif

/****************************************************************************/

/*** Constants ***/

/* minimum allowed density */
#define DENSITY_MIN  1.e-20
#define VDASH_DELTA  1.e-6

/* basic constants */
/*#define PI    M_PI*/
#define PI    3.141592653589793238462643383279502884197169399375

#define KAPPA_H2O 78.5


#define KBOLTZ     1.3807e-23  /* J K-1 */
#define EPSILON_0  8.85419e-12  /* C^2 J^-1 m^-1 */
#define E_CONST    1.60219e-19  /* C */
#define AVOGADRO   6.0220e-23   /* mol^-1 */
#define TRUE  1
#define FALSE 0

#define BINODAL_FLAG  -1325  /* to let initial guess routine know we need to fill X2 */

#define WALL_EL        0 
#define FLUID_EL       1 

#define NCOMP_MAX 5
#define NDIM_MAX  3
#define NWALL_MAX 600 
#define NWALL_MAX_TYPE 20 
#define NPERIODIC_MAX 4

#define NBOND_MAX 4

#define NSTEPS_MAX 10

#define MAX_ROUGH_BLOCK 100
#define FILENAME_LENGTH 300

/* a constat flag to indicate that there is no bond between a pair of segments in 
   a polymer problems */
#define NO_BOND_PAIR -962.0

/* a constant flag for the share_global_vec routine */
#define UNINIT_VEC -200.0

/* constants for first and second derivatives */
#define FIRST  1
#define SECOND 2

/*
 * Types of Nonlinear Solvers available in Tramonto 
 */
#define NEWTON_BUILT_IN       0
#define NEWTON_NOX            1
#define PICARD_BUILT_IN       2
#define PICARD_NOX            3
#define PICNEWTON_BUILT_IN    4
#define PICNEWTON_NOX         5

/*
 * Stencil types refer to the integration schemes needed for different
 * non-local physics. 
 * NSTEN is the total number of stencil types.
 * The following are the current acceptable values for Stencil types.
 */

#define NSTEN        8

#define DELTA_FN_R            0
#define THETA_FN_R            1
#define THETA_PAIRPOT_RCUT    2
#define THETA_CR_RPM_MSA      3
#define THETA_CR_DATA         4
#define THETA_FN_SIG          5
#define DELTA_FN_BOND         6
#define THETA_CR_GENERAL_MSA  7

#define NO_RENORMALIZATION_FLAG -888

#define FROM_LOCA 0
#define FROM_MAIN 1

/*
 * Options for the context flag that is used to assign potential parameters
 */
#define FLUID_FLUID 0
#define WALL_FLUID  1
#define WALL_WALL   2

/*
 * An equation type list is given here to make identification of a particular
 * unknown number straightforward.  This became necessary with the introduction
 * of the WTC polymers where we now have 6 types of equations to fill.
 */
#define NEQ_TYPE       12 
#define NO_UNK        -888

#define DENSITY        0
#define POISSON        1
#define HSRHOBAR       2
#define MF_EQ          3
#define CAVWTC         4
#define BONDWTC        5
#define DIFFUSION      6
#define CMS_FIELD      7
#define WJDC_FIELD     8
#define SCF_CONSTR	   9
#define SCF_FIELD	  10
#define G_CHAIN       11 

/*
#define DENSITY        0
#define CMS_FIELD      1
#define WJDC_FIELD     2
#define SCF_FIELD      3 
#define POISSON        4
#define DIFFUSION      5
#define HSRHOBAR       6
#define CAVWTC         7
#define BONDWTC        8
#define MF_EQ          9
#define G_CHAIN       10 
#define SCF_CONSTR    11
*/


/* Here are some constants needed to make the physics based ordering of the
   matrix an option in the code. */
#define MATRIX_FILL_NODAL 1   /* set to zero for physics based ordering */
#define NODAL_FLAG -999
#define BOX 0
#define LOCAL_N 1
#define GLOBAL 2

/*
 * Zones refer to different quadrature schemes on the mesh.  We allow for 
 * either 1 (all points treated equally) or 3 zones.  In the three zone
 * case the first zone is found within 2 sigma_max of the surfaces.  The
 * second zone is found between 2 and 5 sigma_max of the surfaces, and
 * the third zone is found beyond 5 sigma max of the surfaces 
 */
#define NZONE_MAX  10 

/*
 * These constants identify Surface_Type
 */
#define smooth_planar_wall              0
#define finite_planar_wall              1
#define colloids_cyl_sphere             2
#define atomic_centers                  3
#define point_surface                   4
#define finite_cyl_3D                   5

#define cyl2D_sphere3D_pore             7
#define cyl3D_slit2D_pore               8

/* 
 * These constants identify the functional choices (Type_func).
 */
#define NONE       -1
#define FMT1       0
#define FMT2       1
#define FMT3       2
#define FMT4       3
#define LDA        4
#define GHRM       5
#define GVDWM      6

/*
 *  These constants identify the type of polymer to be studied (Type_poly).
 */
#define NONE        -1
#define CMS          0
#define CMS_SCFT     1
#define WTC          2
#define WJDC         3
#define WJDC2        4 
#define WJDC3        5 
#define SCFT         6	

/*
 *  These constants identify the type of polymer architecture (Type_poly_arch).
 */
#define POLY_ARCH_FILE 0
#define LIN_POLY 1
#define LIN_POLY_SYM 2
#define SET_IN_GUI 3

/*
 * These constants identify attraction functional choices (Type_attr).
 */
#define NONE        -1
#define MFPAIR_RMIN_UMIN      0
#define MFPAIR_SIGMA_ZERO     1
#define MFPAIR_SIGMA_USIGMA   2
#define MFPAIR_SIGTOUMIN_UMIN 3
#define MFPAIR_RCSZERO_ZERO   4

/*
 * These constants set choices for range of inner core of ATT potential (Type_CoreATT_R)
 */
#define ATTCORE_SIGMA       0
#define ATTCORE_UMIN        1
#define ATTCORE_UCSZERO     2
#define ATTCORE_SIGTOUMIN   3

/*
 * These constants set choices for range of inner core of ATT potential (Type_CoreATT_CONST)
 */
#define CORECONST_UCONST    0
#define CORECONST_ZERO      1

/* 
 * These constants identify the functional choices (Type_coul).
 */
#define NONE          -1
#define BARE           0
#define DELTAC_RPM     1 
#define DELTAC_GENERAL 2
#define POLARIZE       3

/*
 * The following are choices for the neutral fluid-fluid interactions (Ipot_ff_n)
 */
#define IDEAL_GAS    0
#define HARD_SPHERE  1
#define LJ12_6       2

/*
 * The following are choices for the charge fluid-fluid interactions (Ipot_ff_c)
 */
#define NO_CHARGE    0
#define COULOMB      1
#define YUKAWA       2

/* The following are choices for Type_vext The VEXT_PAIR_POTENTIAL option uses pair 
   potentials implemented for fluid interactions.  The VEXT_DEFINED option uses 
   external fields defined separately */
#define VEXT_DEFINED  0
#define VEXT_PAIR_POTENTIAL  1

/*
 * The following are choices for pair interacton potentials.  Note that these
 * options can be used to set stencil functions for strict mean field attractions,
 * set external fields based on 3D interaction potentials, and set wall-wall
 * interactions for 3D atomistic surfaces.  They are used by Type_pairPot,
 * Vext_PotentialID, and Type_uwwPot parameters which may be set independently.  
 */
#define PAIR_HARD            -1 
#define PAIR_LJ12_6_CS        0
#define PAIR_COULOMB_CS       1
#define PAIR_COULOMB          2
#define PAIR_YUKAWA_CS        3
#define PAIR_EXP_CS	      4
#define PAIR_SW		      5
#define PAIR_LJandYUKAWA_CS   6
#define PAIR_r12andYUKAWA_CS  7
#define PAIR_r18andYUKAWA_CS  8
#define PAIR_rNandYUKAWA_CS   9


/* options for Type_hsdiam */
#define SIGMA_DIAM          0
#define BH_DIAM             1
#define MANUAL_HS_DIAM         2
/*
 * The following are choices for the Vext_PotentialID parameter.  
 */
#define LJ9_3_CS          0
#define LJ9_3_v2_CS       1
#define LJ9_3_noCS        2
#define LJ9_3_shiftX_CS   3
#define REPULSIVE9_noCS   4
#define EXP_ATT_noCS      5
#define LINEAR_noCS       6 
#define R7_YUKAWA_SUM_CS  7 

/*
 * The following are settings for the physics scaling parameter 
 */
#define AUTOMATIC       1
#define MANUAL_INPUT    2

/*
 * The following are choices for Type_interface which indicates 
 * there the interface requires boundary conditions that differ on
 * two sides of the domain. Note that for a phase interface, density is different while
 * chemical potential is constant (no diffusion).  For a diffusive interface, both density
 * and chemical potential vary, and the diffusion equation must also be solved.
 */
#define UNIFORM_INTERFACE  0
#define DIFFUSIVE_INTERFACE 1
#define PHASE_INTERFACE 2

/*
 *  The following define ways to distribute charge in the system 
 */
#define POINT_CHARGE 0
#define SMEAR_CHARGE 1
#define BACKGROUND   2

/*
 * The following are choices for the initial guess (Iguess_fields)
 */
#define BULK              0
#define CALC_ALL_FIELDS   1
#define CALC_RHOBAR_ONLY  2
#define CALC_RHOBAR_AND_G 3

/*
 * The following are choices for the initial guess (Iguess)
 */
#define CONST_RHO        0 
#define EXP_RHO          1
#define STEP_PROFILE     2
#define CHOP_RHO         3
#define CHOP_RHO_STEP    4
#define LINEAR           5
#define RAND_RHO		 6  //LMH

/*
 * The following are choices for how to handle the density profile
 * on a restart.
 */
#define NORESTART          0
#define RESTART_BASIC      1
#define RESTART_STEP       2
#define RESTART_DENSONLY   3
#define RESTART_FEWERCOMP  4
#define RESTART_1DTOND     5

/*
 * This constant is a flag that allows us to use the fill routines to set up
 * the initial guess.  We need to do this for some of the more complex variables 
 */
#define INIT_GUESS_FLAG  2
#define CALC_RESID_ONLY  3
#define CALC_AND_FILL_RESID_ONLY  4

/*
 * The following are the various fields for continuuation.  There are
 * 3 groups of continuation routines. Group I contains the core Tramonto
 * capabilities.  Group II contains archived extensions and special cases
 * that have been used for specific applications.  Group III should contain
 * local user extensions to be found in dft_plugins_user_continue.c.  
 */
/*
 * This group contains the core continuation capabilties in Tramonto
 */
#define NCONT_MAX          2 /* the maximum number of solutions possible for use with Loca */
#define CONT_MESH          0   /* mesh size */
#define CONT_TEMP          1   /* State Parameters */
#define CONT_RHO_I         2
#define CONT_BETAMU_I      3  /* Vary chemical potential for species I */
#define CONT_EPSW_I        4    /* Wall-Wall Energy Params for wall I */
#define CONT_EPSWF_IJ      5    /* Wall-Fluid Energy Params for IJ term */
#define CONT_EPSFF_IJ      6   /* Fluid-Fluid Energy Params for IJ term */
#define CONT_ELECPARAM_I   7  /* Charged surface params */
#define CONT_ELECPARAM_ALL 8  /* Charged surface params */
#define CONT_SEMIPERM_IJ   9  /* Vext_membrane */
#define CONT_SIGMAFF_IJ    10   /* Fluid-Fluid Interaction Diameter for IJ term */
#define CONT_BETAMU_I_NEW  11 /* Vary chemical potential for species I...holding densities of other species constant */
/*
 * This group contains extensions to the core capabilities that are archived in the repository.
 */
#define CONT_LOG_RHO_I          100
#define CONT_RHO_CONST_RHOTOT58 101
#define CONT_RHO_CONST_XSOLV    102
#define CONT_RHO_ALL		103
#define CONT_EPSW_ALL		104
#define CONT_EPSWF_ALL	        105
#define CONT_EPSWF_SOME 	106
#define CONT_EPSFF_ALL		107
#define CONT_CRFAC              108  /* continuous mixing of two cr files */
/*
 * Any other new contination types that are defined by the user should be set here beginning
 * with the ID number 200
 */

#define PRINT_RHO_0      0

/*
 * These constants identify integration schemes for non-local interactions
 * (AGS) Error checking should probably go in the input file reader, to
 *    make sure that the integration scheme is compatable with the weights.
 */

 /* These choices are for the surface of a sphere*/
#define TETRAHEDRON   1
#define OCTAHEDRON    2
#define CUBE          3
#define ICOSAHEDRON   4
#define DODECAHEDRON  5
#define SEVENTY_TWO   6
#define ONE_D_TWELVE  7
#define MIDPOINT_RULE 8

 /* The following is a choice for the wall type Ipot_wf_n*/
#define VEXT_NONE          0
#define VEXT_HARD        1
#define VEXT_DIST_TO_SURF          2  /* any potential that is a function of only distance from the surface */
#define VEXT_DIST_TO_CENTER        3  /* any potential that is a function of the distance to the center of wall (atom). */
#define VEXT_3D_INTEGRATED      5  /* more proper 3D integration potential for funny geometries */

 /* The following is a choice for the wall type Ipot_ww_n*/
#define NO_WW              0
#define ATOM_CENTERS_WW    1 

 /* The following is a choice for the BC type */
#define IN_WALL             -1
#define IN_BULK              0
#define PERIODIC             1
#define REFLECT              2
#define LAST_NODE            3
#define LAST_NODE_RESTART    4

/* These choices are for types of boundary conditions for
    charged systems                                    */
#define CONST_POTENTIAL  1
#define CONST_CHARGE     2
#define ATOMIC_CHARGE    3

 /* The following is a choice for the load balance flag */
#define LB_LINEAR    0
#define LB_BOX       1
#define LB_WEIGHTS   2

/* The following is a choice for the treatment of dielectric_constants */
#define DIELEC_CONST       0
#define DIELEC_WF          1
#define DIELEC_WF_PORE     2
#define DIELEC_W_DENS      3


/*  The following are choices for a surface element direction from
    a given node.  These are stored in Surf_elem_type[][], and they
    refer to the directions one must move in  3D to get to the
    center of a given surface element whose normal is stored in the
    Surf_normal array !!  In 2D, only the first 4 are relevent, and
    can be read UP,RIGHT,LEFT,DOWN */

/* 2D Example        UP
 *                    |
 *                    |
 *           LEFT --- 0 ---- RIGHT
 *                    |
 *                    |
 *                    DOWN
 */


#define UP_BACK      0
#define RIGHT_BACK   1
#define LEFT_BACK    2
#define DOWN_BACK    3
#define UP_FRONT     4
#define RIGHT_FRONT  5
#define LEFT_FRONT   6
#define DOWN_FRONT   7
#define RIGHT_UP     8
#define LEFT_UP      9
#define RIGHT_DOWN  10 
#define LEFT_DOWN   11 


/* The following are options for 1D steady state transport
   models.  They include: the unit area (OPTION_ONE),
   a cylindrical pore (OPTION_CYL) and a varying diameter
   pore defined by the user (OPTION_VARY) */

#define OPTION_ONE   0
#define OPTION_CYL   1
#define OPTION_VARY  2


/* The following are choices for file output */
#define FILES_BASIC        0
#define FILES_EXTENDED     1 
#define FILES_DEBUG        2
#define FILES_DEBUG_MATRIX 3 

/* The following are choices for screen output */
#define SCREEN_NONE       -1 
#define SCREEN_ERRORS_ONLY  0 
#define SCREEN_BASIC       1
#define SCREEN_DEBUG_RESID 2
#define SCREEN_VERBOSE     3 

/* The following are old choices for output */
#define MINIMAL      0
#define DENSITIES    1 
#define EXTENDED     2
#define VERBOSE      3 
#define NO_SCREEN    4 
#define VERBOSE_MATRIX    5 

/* The followint are choices for the output of density info */
#define SWITCH_NO_STATEOUT -1
#define SWITCH_RHO   0
#define SWITCH_ALLTYPES_ICOMP 1
#define SWITCH_ION   2
#define SWITCH_MU    3
#define SWITCH_ALLTYPES 4
#define SWITCH_BULK_OUTPUT 5
#define SWITCH_BULK_OUTPUT_ALL 6

#define SWITCH_SURFACE_SEP   0

/* The following is a flag for 1D boundary conditions in a 2D or
    3D domain --- only set up for steady state problems now */
#define FLAG_1DBC   -999
#define FLAG_BULK   -888
#define FLAG_PBELEC -777
#define BULK_ZONE    2 
#define PB_ZONE      3 

/* Polymer constants */
#define N_NZCR_MAX   200   /* maximum # of non-zero's in direct correlation fn */
#define NBLOCK_MAX   20 
#define NMER_MAX     200

/* options for reading external field */
#define READ_VEXT_FALSE      0
#define READ_VEXT_TRUE       1
#define READ_VEXT_SUMTWO     2
#define READ_VEXT_STATIC     3

/* options for jacobian coarsening */
#define JAC_RESID_ZONES_SAME    0
#define JAC_ZONE0_FAC2LESSTHANRESID    1
#define JAC_ZONES_FAC2LESSTHANRESID    2
#define JAC_ZONES_ALLMOSTCOARSE        3
#define JAC_ZONES_SECONDMOSTCOARSE     4
#define JAC_ZONES_SETFIXED_ESIZE       5

/* options for how to specify quantity of grafted chains */
#define GRAFT_DENSITY 1
#define GRAFT_NUMBER 2

/****************************************************************************/

/*** Definitions of structures, enumerated types, typdefs, ...  ***/

/* 
 *  Declaration of the Stencil_Struct follows. Info for integrating
 *  non-local terms will be precalculated and stored here. The
 */

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

/* 
 *  Declaration of the RB_Struct, a structure to hold all the Rho_bar
 *  values (scalar and vector) at a given fluid or wall node.
 */

struct RB_Struct {
  double    S0;      /*   1/(4*pi*Ri*Ri) * Delta_fn   */
  double    S1;      /*   1/(4*pi*Ri)    * Delta_fn   */
  double    S2;      /*                    Delta_fn   */
  double    S3;      /*                    Theta_fn   */
  double    V1[NDIM_MAX];      /*  1/(4*pi*Ri) * unit_vec * Delta_Fn   */
  double    V2[NDIM_MAX];      /*                unit_vec * Delta_Fn   */
};

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

/*
 * Declaration of the Trilinos solver stuff follows.
 */
extern void * LinProbMgr_manager;

/* 
 *  Declaration of the Aztec_Struct follows. 
 *  Information for the linear solver is all held here.
 */

/*   DONE WITH THIS STRUCT !! */
/*struct Aztec_Struct {*/
  /* int    options[AZ_OPTIONS_SIZE]; Array used to select solver options.  */
  /* double params[AZ_PARAMS_SIZE];    User selected solver paramters.       */
/*#ifdef DONE_WITH_THESE
  int    proc_config[AZ_PROC_SIZE];* Processor information.                *
  int    *data_org;                * Array to specify data layout          *
  double status[AZ_STATUS_SIZE];   * Information returned from AZ_solve(). *
  int    *update;                  * vector elements updated on this node. *
  int    *external;                * vector elements needed by this node.  *
  int    *update_index;            * ordering of update[] and external[]   *
  int    *extern_index;            * locally on this processor.            *
  int    *bindx;                   * Sparse matrix to be solved is stored  *
  double *val;                     * in these MSR arrays.                  *
  int    N_update;                 * # of unknowns updated on this node    *
  int    nonzeros;                 * # of nonzeros in sparse matrix        *
#endif
};*/

/* 
 *  Declaration of the Loca_Struct follows.
 *  Information for the continuation library is held here.
 */

#define USE_LOCA
struct Loca_Struct {
  int    method;      /* Continuation method                          */
  int    cont_type1;  /* flag specifying the continuation parameter   */
  int    cont_type2;  /* flag specifying the second (free) parameter  */
  int    num_steps;   /* maximum number of continuation steps to take */
  double aggr;        /* step size control parameter                  */
  double step_size;   /* initial continuation step size               */
};

/****************************************************************************/

/*** extern statements for globals ***/


/*extern int    count_zero;
extern int    count_nonzero;
*/

/* Basic Equation info */
extern int Phys2Nunk[NEQ_TYPE];  /* Number of unknowns of a particular equation type */
extern int Phys2Unk_first[NEQ_TYPE]; /* starting unknown number for a given equation type */
extern int Phys2Unk_last[NEQ_TYPE]; /* ending unknown number for a given equation type */
extern int Unk2Phys[3*NCOMP_MAX+2*NMER_MAX+NMER_MAX*NMER_MAX+13]; /* array that gives equation type
                                                                         given an unknown index */
extern int Solver_Unk[3*NCOMP_MAX+2*NMER_MAX+NMER_MAX*NMER_MAX+13]; /* translation to solver ids */

/* Mesh info */


/*************** Global Mesh ********************************/
extern int Print_flag;
extern double **Array_test;

extern int     Ndim;            /* # of spatial dimensions of the problem      */
extern int     Nnodes;          /* # of nodes in the mesh                */
extern int     Nunknowns;       /* # of unknowns in the problem          */
extern int     Nelements;       /* # of elements in the mesh             */
extern int     Elements_plane;  /* # of elements in the x1-x2 plane      */
extern int     Nodes_plane;     /* # of nodes in the x1-x2 plane          */
extern int Nodes_x[NDIM_MAX];   /* Array[Ndim]: # nodes in each dimension */
extern int Elements_x[NDIM_MAX];/* Array[Ndim]: # elements in each dim    */
extern int     Max_sten_length[3];  /* The number of nodes in the longest stencil */

/*************** Reference Variables ***********************/
extern int Open_GUI;
extern double Length_ref;
extern double Density_ref;
extern double Dielec_ref;
extern double Potential_ref;
extern double VEXT_MAX;

/*************** Extended Local Mesh ************************/
extern int     Nnodes_box;          /* # of nodes in the extended local mesh   */
extern int     Nunknowns_box;       /* # of unknowns in the extended local mesh */
extern int     Nelements_box;       /* # of elements in the extended local mesh */
extern int     Elements_plane_box;  /* # of el in x1-x2 plane  of extended local mesh  */
extern int     Nodes_plane_box;     /* # of nodes in x1-x2 plane of extended local mesh*/
extern int Nodes_x_box[NDIM_MAX];   /* Array[Ndim]: # nodes in idim on extended local mesh */
extern int Elements_x_box[NDIM_MAX];/* Array[Ndim]: # elements in idim on extended local mesh */
extern int     Min_IJK_box[3]; /* The minimum IJK of the extended internals on this procesor */
extern int     Max_IJK_box[3];  /* The maximum IJK of the extended internals on this procesor */
extern int     Pflag[3];   /* PERIODIC flag array[Ndim]: TRUE Nodes_x_box=Nodes_x */

/************** Local Mesh **********************************/
extern int     Min_IJK[3];      /* The minimum IJK of the internals on this procesor */
extern int     Max_IJK[3];      /* The maximum IJK of the internals on this procesor */

/************** Communications arrays for Gather vectors ******/
extern int *Comm_node_proc;   /* array on proc 0 of nodes per processors */
extern int *Comm_unk_proc;   /* array on proc 0 of unknowns per processors */
extern int *Comm_offset_node; /* array on proc 0 of offsets of nodes  MPI_Gatherv*/
extern int *Comm_offset_unk;  /* array on proc 0 of offsets of unknowns MPI_Gatherv*/

/************** Other **************************************/
extern double  Vol_el;          /* Volume of one element of our regular mesh        */
extern double  Area_surf_el[3]; /*Area of surface element with normal in idim direction*/
extern double Vol_in_surfs[NCOMP_MAX];  /* volume in all of the surfaces for each list */
extern  int     Nlists_HW;       /* Number of lists needed if hard walls (for mixtures)*/
extern int     Nel_wall_count;  /* Number of elements in the 0.5sigma of wall:rhobars*/
extern int    **Nelems_per_wall;   /* Number of elements in a given [iwall] wall        */
extern int     Nnodes_per_proc;  /* Number of nodes owned by this processor         */
extern int     Nunk_int_and_ext; /* Number of unknownsneeded on this processor      */
extern int     *B2L_unknowns;     /* Box to Local array for all unknowns */
extern int     *B2G_node;         /* Box to global array for all box nodes */
extern int     *B2G_unk;          /* Box to global array for all box unknowns */
extern int     *L2B_node;         /* Local to box array for all local nodes */
extern int     *B2L_node;         /* Box to local array for all local nodes */
extern int     *L2G_node;         /* Local to global coordinates */
extern int     Nunk_per_node;   /* Number of unknowns per node (usually Ncomp       */
extern int     Nmf_eqns;        /* Number of mean field - attraction equations */
extern int     Nrho_bar;        /* Number of rhobar equations per node */
extern int     Nrho_bar_s;      /* Number of scalar rhobar equations per node */
extern int     Npoisson;        /* Number of rhobar equations per node */
extern int     Ndiffusion;        /* Number of rhobar equations per node */
extern int     Ndensity_unk;         /* Number of unknowns for the Euler-Lagrange equation */
extern int     Ntype_unk;       /* Number of equations defining segment types for polymer TC cases */
extern int     Nrho_bar_cavity; /* Number of nonlocal densities for the cavity function - WTC polymers */
extern int     Nrho_bar_bond; /* Number of nonlocal densities for the bond functionals - WTC polymers */
extern int     Nnodes_per_el_V;  /* Number of nodes per volume element              */
extern int     Nnodes_per_el_S;  /* Number of nodes per surface element            */
extern int     Plane_new_nodes; /* Indicates in which plane (xy,yz,xz)nodes are added*/
extern int     Pos_new_nodes;   /* Indicates where nodes are added                  */
extern double  Size_x[NDIM_MAX];    /*Array of the size of the domain in each dim. */
extern double  Esize_x[NDIM_MAX];   /*Array of the size of an element in each dim. */
extern int     Lmesh_refine;       /*Switch for auto mesh refinement               */
extern int     Type_bc[NDIM_MAX][2];/*Array of boundary conditions in each dim.    */
extern int     Non_unique_G2B[4];   /* Flag indicating which box dimensions       
                                wrap around completely, resulting in a non-unique G2B mapping */
extern int   **Nodes_2_boundary_wall;  /*Array[Nlists_HW][Nnodes] -1 if nod b.n. else b.n.  */
extern int   **Wall_elems;     /*Array[Nlists_HW][Nelements] TRUE for wall elements */
extern int   ****Touch_domain_boundary; /*Array[Nwall][Nlists_HW][[Ndim][2] =0 if surface hits left
                                  boundary = 1 if hits right domain boundary else -1 */
extern int   **Xtest_reflect_TF; /* Array[Nwall][Ndim] for reflections/wall-wall boundaries */
/* these are used to set up the Nodes_2_boundary_wall array then discarded*/
extern int   Nnodes_wall_box; /* Count number of nodes in box that touch a wall */
extern int  *Nodes_wall_box; /* Array to store which nodes touch a wall */
extern int  *Nwall_touch_node; /* Array to store number of walls touching a given node */
extern int  **Wall_touch_node; /*Array to store which walls touch a given node */
extern int  **List_wall_node; /*Array to store which walls touch a given node */
extern int  *Index_wall_nodes; /* ArraY to store indexing in these mesh arrays */

extern int First_time; /* for MSR preprocessing */

extern int     Nzone;          /* Number of diff. quadrature zones on the mesh      */
extern int    *Nodes_to_zone;   /* Array[Nnodes] of quadrature zones */
extern double  Rmax_zone[5];    /* Array distances from surfaces in quadrature zones */
extern int     Mesh_coarsening;  /* Flag indicating whether mesh coarsening is on */
extern int    *Mesh_coarsen_flag;/* Flag (Nnodes) telling how much coarsening for
                             a given node, or negative values telling which*/
extern int     Nnodes_coarse_loc; /* Number of coarse nodes local to a processor */
extern int   *List_coarse_nodes; /* List of coarse nodes local to a processor */
extern int     Coarser_jac;     /* Flag to switch on coarser jacobian than residual */
extern double  Jac_grid;     /* Flag to switch on coarser jacobian than residual */
extern int     Lcut_jac;  /* Logical to indicate if Jacobian stencils will be cut off */
extern double  Jac_threshold; /* Threshold level for Jacobian stencils ... max/Jac_threshold */
extern int    Constant_row_flag[NEQ_TYPE]; /* A flag to turn off calculation of constant coefficients after first fill */

/* Continuation info */
extern int     Nodes_old;         /* # of nodes in the mesh of the previous run */
extern int Nodes_x_old[NDIM_MAX];/* Array[Ndim]: # nodes in idim of previous run */
extern double  *X_old;           /* Array of unknowns from previous run */
extern double  *X2_old;           /* Array of unknowns from previous run */

extern int     Print_force_type;  /* flag for printing of force */
extern int     Print_rho_type;  /* flag for file type for printing of densities */
extern int     Print_rho_switch; /* flag for printing densities -- format */
extern int     Lprint_gofr; /* flag for printing radial distribution functions */
extern int     Lprint_scaleFacWJDC;  /* flag for printing ScaleFac array */
extern int     Lprint_pmf; /* flag for printing radial distribution functions */
extern int     Print_mesh_switch; /* flag for printing densities -- format */
extern int     Lper_area;  /*logical for per unit are outputs of params */ 
extern int     Lcount_reflect;  /*logical for per unit are outputs of params */ 
extern int     Nruns;           /* Number of runs to perform (varying the mesh)     */
extern double  Del_1[NWALL_MAX_TYPE];    /*Stepping parameter for field #1  */
extern double  Rho_max;         /* max rho when using an old solution for mesh contin*/
extern int     Imain_loop;    /* Couter on the number of times through the program  */


/* Surface Physics info */
extern int     Nwall;           /* Number of surfaces in the calculation            */
extern int     Nwall_type;       /* Number of surface types in the problem */
extern int     WallType[NWALL_MAX]; /* array containing type number for each surface */
extern int     Nlink;           /* Number of macro-surfaces */
extern int     Link[NWALL_MAX]; /* Index iwall to a particular macro-surface */
extern int     *Nwall_this_link; /* Number of walls linked together in a particular macro-surface */
extern int     **Link_list; /* List of walls linked in a particular macro-surface */
extern int     Surface_type[NWALL_MAX_TYPE];    /* Type of surfaces of interest                     */
extern int     Orientation[NWALL_MAX_TYPE];  /* Orientation of planar/bumpy infinite walls*/
extern double  WallParam[NWALL_MAX_TYPE];/* Array[Nwall] of a characteristic wall parameter*/
extern double  WallParam_2[NWALL_MAX_TYPE];/* Array[Nwall] of a characteristic wall parameter*/
extern double  WallParam_3[NWALL_MAX_TYPE];/* Array[Nwall] of a characteristic wall parameter*/
extern int     Lapply_offset[3];/* Array of logicals to control how the offsets are applied to various WallParams*/
extern int     Lrough_surf[NWALL_MAX_TYPE]; /*Logical for rough surfaces */
extern double  Rough_precalc[NWALL_MAX_TYPE][MAX_ROUGH_BLOCK][MAX_ROUGH_BLOCK];
extern double  Rough_length[NWALL_MAX_TYPE];
extern double  Rough_param_max[NWALL_MAX_TYPE];
extern int     read_rough; /* a way to indicate whether surface roughness params are read in */
extern int     read_periodic; /* a way to indicate whether periodic surfaces are being used */
extern int     read_wedge; /* a way to indicate whether wedge cutouts are being used */
extern int     read_linear; /* a way to indicate whether linear surface modifications are being used */
extern int     Lwedge_cutout[NWALL_MAX];    /* TRUE or FALSE for applying wedge cutout to surface*/
extern double  Angle_wedge_start[NWALL_MAX];  /* start angle for wedge cutout */
extern double  Angle_wedge_end[NWALL_MAX];     /* end angle for wedge cutout */
extern  int    Lperiodic_overlay[NWALL_MAX_TYPE];    /* TRUE or FALSE for periodic function added to surface */
extern  int    Nperiodic_overlay[NWALL_MAX_TYPE];     /* The number of periodic functions to apply */
extern  int    OrientationPeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];     /* The orientation of periodic functions to apply */
extern  double AmplitudePeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];     /* The amplitude of periodic functions to apply */
extern  double WavelengthPeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];    /* The period of periodic functions to apply */
extern  double OriginPeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];     /* The origin of periodic functions to apply */
extern  int    Llinear_overlay[NWALL_MAX_TYPE];    /* TRUE or FALSE for linear function added to surface */
extern  int    Nlinear_overlay[NWALL_MAX_TYPE];     /* The number of linear functions to apply */
extern  int    OrientationLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];     /* The orientation of linear functions to apply */
extern  double SlopeLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];     /* The slope of linear functions to apply */
extern  double OriginLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];     /* The origin of linear functions to apply */
extern  double EndpointLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];     /* The endpoint of linear functions to apply */
extern int     Lauto_center; /* Logical to automatically center the surfaces in the domain */
extern int     Lauto_size; /* Logical to automatically size the domain */
extern int     Lrandom_walls; /* Logical to turn on random wall placement */
extern double  WallPos[NDIM_MAX][NWALL_MAX]; /* Array of the centers of the surfaces*/
extern double  **WallPos_Images; /* Array of the centers of the surfaces including periodic and reflected images*/
extern int     *WallType_Images; /* Array of the types of the surfaces including periodic and reflected images*/
extern int     *RealWall_Images; /* Array of the real wall in the domain associated with a particular images*/
extern int     *Image_IDCheck; /* Array of the real wall in the domain associated with a particular images*/
extern int     Nwall_Images; /* Number of surfaces including all images*/

/* Fluid Physics info */
extern int     Ncomp;           /* Number of components in the current problem      */
extern double  Temp;            /* Reduced Temperature for LJ interactions          */
extern double  Temp_elec;       /* Reduced Temperature for Poisson's equation       */
extern double  charge_fluid;    /* the charge in the fluid ... post processing */
extern  int     Type_dielec;      /* choose how to handle dielectric constants in system */
extern double  Sigma_Angstroms_plasma;     /* the particle diameter in Angstroms sigma to be used in computing the plasma parameter. */
extern double  Temp_K_plasma;     /* the temperature in K to be used in computing the plasma parameter.*/
extern double  DielecConst_plasma;     /* the dielectric constant to be used in computing the plasma parameter */
extern double  Dielec_bulk;     /* the dielectric constant in the bulk fluid */
extern double  Dielec_pore;     /* the dielectric constant in the "pore" fluid */
extern double  Dielec_X;        /* distance from a surface defines the "pore" fluid */
extern double *Dielec_wall;     /* the dielectric constant as function of wall type */
extern double *Dielec;     /* dielectric constant as function of ielement_box */
extern int     Lpolarize[NCOMP_MAX]; /* logical for if a species is polarizeable */
extern double  Pol[NCOMP_MAX];  /* bulk polarizeability for each species */
extern double  Energy;   /* Surface free energy to return to Towhee */


extern double  Betap;           /* Pressure in units of kT sigma_ff[1]^3           */
extern double Betap_LBB;       /* Pressure calculated for LBB of domain */
extern double Betap_RTF;       /* Pressure calculated for RTF of domain */
extern double  Betap_id;       /* Ideal gas Presseure in units of kT sigma_ff[1]^3   */
extern double  Betap_att;      /* Attractive Presseure in units of kT sigma_ff[1]^3   */
extern double  P_over_po;
extern int     L_isotherm; /* Logical for isotherm vs. force per distance data */
extern double  Rho_b[NCOMP_MAX];   /* Array[Ncomp] of component bulk densities      */
extern double  Rho_t;
extern double  Rhobar_b[10]; /* Array[Nrho_bar] of bulk rhobars      */
extern double  Rhobar_b_LBB[10]; /* Array[Nrho_bar] of bulk rhobars LBB  */
extern double  Rhobar_b_RTF[10]; /* Array[Nrho_bar] of bulk rhobars RTF  */
extern double  Dphi_Drhobar_b[10]; /* Array[Nrho_bar] of bulk energy derivs w/r/to rhobars      */
extern double  Dphi_Drhobar_LBB[10]; /* Array[Nrho_bar] of bulk energy derivs w/r/to rhobars LBB  */
extern double  Dphi_Drhobar_RTF[10]; /* Array[Nrho_bar] of bulk energy derivs w/r/to rhobars RTF  */
extern double  Rho_seg_b[NMER_MAX]; /* array of bulk segment densities */
extern double  Rho_seg_LBB[NMER_MAX];
extern double  Rho_seg_RTF[NMER_MAX];
extern double Field_WJDC_b[NMER_MAX];
extern double Field_WJDC_LBB[NMER_MAX];
extern double Field_WJDC_RTF[NMER_MAX];
extern double Field_CMS_b[NMER_MAX];
extern double Field_CMS_LBB[NMER_MAX];
extern double Field_CMS_RTF[NMER_MAX];
extern double Scale_fac_WJDC[NCOMP_MAX][NCOMP_MAX];
extern double G_WJDC_b[NMER_MAX*NBOND_MAX];
extern double G_WJDC_LBB[NMER_MAX*NBOND_MAX];
extern double G_WJDC_RTF[NMER_MAX*NBOND_MAX];
extern double G_CMS_b[NMER_MAX*NBOND_MAX];
extern double G_CMS_LBB[NMER_MAX*NBOND_MAX];
extern double G_CMS_RTF[NMER_MAX*NBOND_MAX];
extern double  *Rhobar3_old;   /* Array[Nnodes_box] of old values of rhobar 3*/
extern double Xi_cav_b[4]; /* Array of bulk rhobars for cavity functions of WTC polymer functionals */
extern double Xi_cav_LBB[4]; /* Array of bulk rhobars for cavity functions of WTC polymer functionals */
extern double Xi_cav_RTF[4]; /* Array of bulk rhobars for cavity functions of WTC polymer functionals */
extern double BondWTC_b[NMER_MAX*NMER_MAX]; /*Array of bulk rhobars for bonds in WTC functionals*/
extern double BondWTC_LBB[NMER_MAX*NMER_MAX]; /*Array of bulk rhobars for bonds in WTC functionals*/
extern double BondWTC_RTF[NMER_MAX*NMER_MAX]; /*Array of bulk rhobars for bonds in WTC functionals*/
extern double  Rho_coex[2];   /* Liquid and Vapor Coexisting Densities         */
extern double  Betamu_hs_ex[NCOMP_MAX];/* Array of excess hardsphere chemical potentials*/
extern double  Betamu[NCOMP_MAX];   /* Array[Ncomp] of chemical potentials*/
extern double  Betamu_id[NCOMP_MAX];   /* Array[Ncomp] of ideal gas chemical potentials*/
extern double  Betamu_wtc[NMER_MAX];
extern double  Betamu_wtc_LBB[NMER_MAX];
extern double  Betamu_wtc_RTF[NMER_MAX];
extern double  Betamu_chain[NMER_MAX];
extern double  Betamu_chain_LBB[NMER_MAX];
extern double  Betamu_chain_RTF[NMER_MAX];
extern double  Betamu_ex_bondTC[NCOMP_MAX][NMER_MAX*NMER_MAX];/* Array of excess segment chemical potentials - WTC poolymer*/
extern double  Betamu_seg[NMER_MAX];/* Array of excess segment chemical potentials - WTC poolymer*/
extern double  Betamu_seg_LBB[NMER_MAX];/* Array of excess segment chemical potentials - WTC poolymer*/
extern double  Betamu_seg_RTF[NMER_MAX];/* Array of excess segment chemical potentials - WTC poolymer*/
extern int     Ipot_ff_n;    /* Potential Type for neutral part of f-f interactions */
extern int     Ipot_wf_n[NWALL_MAX_TYPE];    /* Potential Type for neutral part of w-f interactions */
extern int     Type_pairPot;  /* Interaction potential to use for strict mean field DFT calculations*/
extern int     Type_hsdiam;  /* How to calculate Hard-Sphere diamters - use Sigma or Barker-Henderson approach */
extern int     Type_vext[NWALL_MAX_TYPE];  /* External field type for a given surface type - either pair potential, u(r) or defined vext*/
extern int     Vext_PotentialID[NWALL_MAX_TYPE];  /* ID for external field for each surface type in the problem */
extern int     Type_uwwPot;  /* potential to use for computation of wall-wall interactions.  Used in 3D-atomic surface calculations */
extern int     Ipot_ww_n[NWALL_MAX_TYPE][NWALL_MAX_TYPE];    /* Potential Type for neutral part of w-f interactions */
extern int     Ipot_ff_c;    /* Potential Type for charged part of f-f interactions */
extern int     Ipot_wf_c;    /* Potential Type for charged part of w-f interactions */
extern int     Lhard_surf;   /* Logical indicating if the surfaces have hard cores */
extern int     Lvext_dash;   /* Logical indicating if the Vext_dash array should be set up */
extern int     Iguess;        /* Type of initial guess */
extern int     Iguess_fields;        /* Type of initial guess */
extern double  random_rho;			/*Amount of randomness to add to rho for random initial guess LMH*/
extern int     Nsteps;         /* Number of steps for a step profile initial guess */
extern int     Orientation_step[NSTEPS_MAX]; /* orientation of the step profile */
extern double  Xstart_step[NSTEPS_MAX];  /* start position array for the step profile */
extern double  Xend_step[NSTEPS_MAX];  /* end position array for the step profile */
extern double  Rho_step[NCOMP_MAX][NSTEPS_MAX];  /* density array for a step profile */
extern int     Lbinodal;        /* Logical TF for binodal calculation */
extern double  Thickness;    /* Thickness parameter for doing wetting studies */
extern int     Mix_type;  /* Choice of Mixing Rules */
extern double  Mass[NCOMP_MAX];           /* Array of the mass of each specie*/
extern double  Sigma_ff[NCOMP_MAX][NCOMP_MAX];/* Array of f-f interaction diameters */
extern double  Bond_ff[NCOMP_MAX][NCOMP_MAX];/* Array of f-f bond lengths for polymers */
extern double  Fac_overlap[NCOMP_MAX][NCOMP_MAX];/* Array of f-f bond lengths for polymers */
extern double  Fac_overlap_hs[NCOMP_MAX];/* Array of f-f bond lengths for polymers */
extern double  Eps_ff[NCOMP_MAX][NCOMP_MAX];  /* Array of f-f interaction energies  */
extern double  Cut_ff[NCOMP_MAX][NCOMP_MAX];  /* Array of f-f cutoff distances      */
extern double  Npow_ff[NCOMP_MAX][NCOMP_MAX]; /* array of N for r^N potentials */
extern double  Rmin_ff[NCOMP_MAX][NCOMP_MAX];  /* Array of f-f distances to the minimum of a pair potential.      */
extern double  Rzero_ff[NCOMP_MAX][NCOMP_MAX];  /* Array of f-f distances to the location where 
                                                   the cut and shifted pair potential is zero.      */
extern double  Charge_f[NCOMP_MAX];           /* Array of the valence of each specie*/
extern double  Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];/* Array of w-f interaction diameters */
extern double  Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];  /* Array of w-f interaction energies  */
extern double  Cut_wf[NCOMP_MAX][NWALL_MAX_TYPE];  /* Array of w-f cutoff distances      */
extern double  Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];/* Array of w-w interaction diameters */
extern double  Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];  /* Array of w-w interaction energies  */
extern double  Cut_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];  /* Array of w-w cutoff distances      */
extern double  EpsYukawa_ff[NCOMP_MAX][NCOMP_MAX]; /* Yukawa prefactor for fluid-fluid interactions */
extern double  EpsYukawa_wf[NCOMP_MAX][NWALL_MAX_TYPE]; /* Yukawa prefactor for wallfluid interactions */
extern double  EpsYukawa_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE]; /* Yukawa prefactor for wall-wall interactions */
extern double  YukawaK_ff[NCOMP_MAX][NCOMP_MAX]; /* Yukawa decay constant on fluid-fluid interactions */
extern double  YukawaK_wf[NCOMP_MAX][NWALL_MAX_TYPE]; /* Yukawa decay constant on wall-fluid interactions */
extern double  YukawaK_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE]; /* Yukawa decay constant on wall-wall interactions */

extern int     **Lsemiperm; /*Array of logicals for semi-permeable surfaces */
extern double  **Vext_membrane; /*Array potentials for semi-perm surfaces */
extern double  **Vext_set;      /*Array of maximum set points for ext potentials */
extern double **Vext;        /* External field array [Nnodes_local][Ncomp]           */
extern double **Vext_static; /* Static part of external field [Nnodes_local][Ncomp] */
extern double **Vext_coul;        /* Coulomb External field array [Nnodes]           */
extern double *Vext_old;        /* for post processing: ext field array           */
extern double ***Vext_dash;  /* derivative of external field [Nnodes][Ncomp][Nwall]*/
extern double **Uww;        /* Wall-Wall interactions [Nwall-1][Nwall-1]           */
extern double **Uww_link;        /* Wall-Wall interactions [Nlink-1][Nlink-1]           */
extern double **X_wall;  /* Distance from inode to iwall [Nnodes][Nwall]    */
extern double ***Xwall_delUP;  /* Distance from inode+delta to iwall [Nnodes][Nwall]  ... used to compute Vdash  */
extern double ***Xwall_delDOWN;  /* Distance from inode-delta to iwall [Nnodes][Nwall] ... used to compute Vdash   */
extern int    **Zero_density_TF; /* array [Nnodes][icomp] indicates where VEXT_MAX */
extern double  Betamu_att[NCOMP_MAX];   /* sum over jcomp of Van der waals constant a(icomp,jcomp)*/
extern double  Avdw[NCOMP_MAX][NCOMP_MAX];    /*  Van der waals constant a(icomp,jcomp)*/

extern double  Rho_w[NWALL_MAX_TYPE];    /* Array[Nwall] of w-w interaction energies     */
extern double  Elec_param_w[NWALL_MAX]; /* Array: surf charge(potential) per area  */
extern int     Type_bc_elec[NWALL_MAX_TYPE];/* Array of surface b.c.'s for charged systems    */
extern int     Nlocal_charge; /*Number of localized surface charges in the system */
extern double  **Charge_x;    /*Position of a given local charge [Nlocal_charge][Ndim]*/
extern double  *Charge_Diam;  /*Diameter of a given local charge [Nlocal_charge]*/
extern double  *Charge;  /*Value of the local charge [Nlocal_charge]*/
extern int     Charge_type_atoms; /* Type of charge distribution on the atoms */
extern int     Charge_type_local; /* Type of charge distribution on the added local charge */
extern double  *Deltac_b;   /* Array [icomp] of electrostatic correlations in bulk*/
extern double  X_MSA[NCOMP_MAX]; /* Array needed for general MSA electrostatics Oleksy and Hanson */
extern double  Gamma_MSA; /* Parameter needed for general MSA electrostatics Oleksy and Hanson */
extern double  N_MSA[NCOMP_MAX]; /* Array needed for general MSA electrostatics Oleksy and Hanson */
extern double  MSAgen_term1[NCOMP_MAX][NCOMP_MAX]; /* Array needed for general MSA electrosattics Oleksy and Hanson */
extern double  MSAgen_term2[NCOMP_MAX][NCOMP_MAX]; /* Array needed for general MSA electrosattics Oleksy and Hanson */
extern double  MSAgen_term3[NCOMP_MAX][NCOMP_MAX]; /* Array needed for general MSA electrosattics Oleksy and Hanson */
extern double  MSAgen_term4[NCOMP_MAX][NCOMP_MAX]; /* Array needed for general MSA electrosattics Oleksy and Hanson */
extern double  **Charge_w_sum_els; /*Array[Nnodes_b][Ndim] of surface charge per area*/
extern double  *Charge_vol_els; /*Array[Nelemts_box] of volume charge per element */
extern int     Vol_charge_flag; /* Flag for volumetric charges */
extern int     Surf_charge_flag; /* Flag for volumetric charges */
extern int     **Nelems_S; /* Array[Nlists_HW][Nnode_per_proc] of # surf elems the b.node touches*/
extern int     ***Surf_normal; /*Array[Nlists][Nnodes_per_proc][Nelems_S] of unit normal vectors */
extern int     ***Surf_elem_to_wall; /*Array of wall to which a given surface element belongs...
                                 a given node may belong to more than one wall !! */
extern int     **Surf_elem_type; /*Array[Nnodes_per_proc][Nelems_S] of surf elem type */
extern double  ***S_area; /*Array[Nlists][Nwall][Ndim] of total surf. area in idim on iwall*/ 
extern double  **S_area_tot; /*Array[Nlists][Nwall] total surf. area on iwall*/ 
extern int     **Nwall_owners; /*Array[Nilists][el_box] of number of walls 
                         (including images) that own a given element */
extern int     ***Wall_owners; /*Array[ilist][iel_box][Nwall_owners] that stores
                          all of the wall owners of a given element */ 

/* surface arrays for cases where global surface information is required */
extern int *NodesS_global;   /* total number of surface nodes in the problem */
extern int **NodesS_GID_global;   /* Global IDs for every surface nodes in the problem */
extern int **NelemsS_global;  /* store the number of surface elements touched by each surface node */
extern int ***Surf_normal_global;  /* store the surface normals for each surface element touched by a surface node */
extern int ***Surf_elem_to_wall_global;  /* store the wall ID associated with each surface element touched by each surface node */
extern int Nnodes_box_extra;  /* Nnodes_box augmented by surface nodes */
extern int *B2G_node_extra; /* array to index Box to global for all nodes including the extra surface nodes */
extern int **S2B_node; /* array to translate surface nodes to box nodes for efficient Gsum calculations*/

/* Steady State Solutions Info */
extern int    Type_interface;          /*Select type if interfacial problem to study*/
extern int    Lconstrain_interface;   /*Logical to control interface constraint*/
extern int    LBulk;          /*True-False Indicates a bulk run - changes output*/
extern int    Linear_transport;       /*True-False Steady State or Equilibrium Run*/
extern double Velocity;               /*Constant Convective Velocity over Diffusion coefficient*/
extern int    Grad_dim;               /*direction where gradient is implemented*/
extern int    Dim_1Dbc;               /*direction where we expect to have a 1D solution near the boundaries*/
extern int    L1D_bc;               /*logical for 1D boundary condition in Grad_dim direction */
extern double X_1D_bc;               /*distance where 1D boundary is applied */
extern double X_const_mu;             /*length where constant mu applies @ edges of domain*/
extern double Rho_b_LBB[NCOMP_MAX];   /*Rho_b boundary condition left-bottom-back*/
extern double Rho_b_RTF[NCOMP_MAX];   /*Rho_b boundary condition right-top-front*/
extern double Elec_pot_LBB;           /*Electric potential boundary condition LBB*/
extern double Elec_pot_RTF;           /*Electric potential boundary condition RTF */
extern int    Flag_mV_elecpot;         /* TF logical for units entry for electrostatic potential boundary conditions */
extern double Betamu_LBB[NCOMP_MAX];  /*Chemical Potential Boundary Condition LBB */
extern double Betamu_RTF[NCOMP_MAX];  /*Chemical Potential Boudary Condition RTF*/
extern double D_coef[NCOMP_MAX];  /*Diffusion Coefficients for ion species */
extern double *Pore_rad_L_IC;    /* array of left  Radii of ion chan pore segments (1D) */
extern double *Pore_rad_R_IC;    /* array of right Radii of ion chan pore segments(1D) */
extern double *Lseg_IC;          /* array of length of ion chan pore segments (1D) */
extern int    Nseg_IC;           /* number of pore segments in a given ion channel */
extern double    *Area_IC;      /* 1D ion channel area per node (box units)*/
extern int   Geom_flag;    /* geometry flag for ion chan. see OPTION_ definitions*/

/* OUTPUT INTEGRAL PARAMETERS */
extern int    **Nel_hit;      /* number of elements hit by a given node in a given list */
extern int    **Nel_hit2;     /* same as prev. for a bulk fluid */
extern int    List[2];       /* which list numbers we care about for integrals*/
extern int    Imax;          /* how many lists are relevent to the case at hand 1 or 2 */
extern double Area;
extern double Fac_vol;
extern double Fac_area;

/* SOME CONSTANTS */
extern double Inv_4pi;               /* Precalculated value of 1/(4*pi)                    */
extern double Inv_rad[NCOMP_MAX];    /* Precalculated inverse of component radius          */
extern double Inv_4pir[NCOMP_MAX];   /* Precalculated inverse of component's 4 pi Radius   */
extern double Inv_4pirsq[NCOMP_MAX]; /* Precalculated inverse of component's 4 pi Radius^2 */

extern int     Type_func;    /* Type of functional for the calculation              */
extern int     Type_attr;    /* Type for handling attractions                       */
extern int     Type_CoreATT_R;     /* Type for range of constant core region for attractions */
extern int     Type_CoreATT_CONST;  /* Type for value of constant in constant core region for attractions */
extern int     Type_coul;    /* Type for handling coulomb interactions              */
extern int     Type_poly;    /* Type for handling polymers                          */
extern int     Type_poly_arch;    /* Type of polymer architecture                        */
extern int     Lseg_densities; /* Logical to indicate that segement (rather than component) densities are treated in the code */
extern int     L_HSperturbation; /* Logical to indicate whether the run is base on perturbation of hard spheres */
extern int     LDeBroglie; /* logical to turn on the DeBroglie wavelength contribution to the free energy functional */

/* Hard core type */
extern double      HS_diam[NCOMP_MAX];  /* Hard sphere diameters for the calculation */

/* Startup Info */
extern int     Restart;     /* Logical that switches between new prof & restart file*/
extern int     Restart_field[NEQ_TYPE];
extern int     Nmissing_densities; /* special restart case where only partial densities are in restart file */
extern int     Restart_Vext;     /* Logical that defines reading of external field*/
extern char    *Vext_filename;       /* pointer to the vext filename */
extern char    *Vext_filename2;       /* pointer to the vext2 filename */
extern char    vext_file_array[FILENAME_LENGTH];       /* file name that contains external field to read in */
extern char    vext_file2_array[FILENAME_LENGTH];       /* a second file name that contains another part of the external field to read in */
extern int     Iwrite;       /* Do we want a complete or modified set of output data*/
extern int     Iwrite_screen;       /* Do we want a complete or modified set of output data*/
extern int     Iwrite_files;       /* Do we want a complete or modified set of output data*/

/* Parallel Info, Aztec info */
extern int     Num_Proc; /* The total number of processors used in this calculation */
extern int     Proc;     /* The unique  processor number (from 0 to Num_Proc-1)     */
extern struct  Aztec_Struct Aztec; /* Structure to hold all the Aztec info          */
extern int     Load_Bal_Flag; /* Flag specifying type of laod balancing to do       */
extern int L_Schur; /* Switch to turn on Schur solvers */
extern void * ParameterList_list; /* Parameterlist to hold Aztec options and params info */


/* Nonlinear Solver info */
extern int NL_Solver;    /* select type of nonliear solver */
extern int Max_NL_iter;    /* Maximum # of Newton iterations (10 - 30)          */
extern int Physics_scaling; /* do physical scaling of nonlinear problems */
extern int ATTInA22Block; /* Logical for location of dense attractions.  1=TRUE=A22block; 0=FALSE=A12block */
extern int Analyt_WJDC_Jac; /* Logical for handling of WJDC jacobians - 0=FALSE=approximate jacobian; 1=TRUE=analytic */
extern double NL_abs_tol,NL_rel_tol; /* Convergence tolerances (update_soln)*/
extern double NL_abs_tol_picard,NL_rel_tol_picard; /* Convergence tolerances (update_soln) --- may be different than newton tolerances*/
extern double NL_update_scalingParam; /* Minimum fraction to update solution to slow down
                           Newton's method */

/* Timers */
extern double Time_linsolver_first;
extern double Time_linsolver_av;
extern double Time_manager_first;
extern double Time_manager_av;
extern double Time_fill_first;
extern double Time_fill_av;
extern double Time_NLSolve;
extern double Time_InitGuess;

/* Linear Solver info */
extern int Az_solver;
extern int Az_kspace;
extern int Az_scaling;
extern int Az_preconditioner;
extern double Az_ilut_fill_param;
extern double Az_tolerance;
extern int    Max_gmres_iter;


extern struct  Loca_Struct Loca; /* Information for continuation library */
/*
 * The global variable Stencil is a 3D array of Stencil_Struct
 * (for hard spheres of size [Nsten][Ncomp][Ncomp] ). This variable is
 * defined in dft_stencil.h, and the extern statement is found here.
 */

extern struct Stencil_Struct ***Stencil;
extern struct SurfaceGeom_Struct *SGeom;

extern int MPsten_Npts_R[NZONE_MAX];  /* # radial gauss pts. in MIDPOINT rule */
extern int MPsten_Npts_arc[NZONE_MAX]; /* # theta gauss pts. in MIDPOINT rule */
extern int MPsten_Npts_phi[NZONE_MAX]; /* # phi gauss pts. in MIDPOINT rule */
extern int Sten_Type[NSTEN]; /* on/off Flag for stencil types  */
extern int Sten_Choice_S[NSTEN][NZONE_MAX]; /* Quadrature type fore each stencil   */
extern int Sten_Choice_R[NSTEN][NZONE_MAX]; /* # Radial Gauss points for THETA_FNCs */

/* Polymer variables */
extern double Deltar_cr;
extern double ***Rism_cr;
extern double Crfac;
extern double Cr_rad[NCOMP_MAX][NCOMP_MAX];
extern double Cr_rad_hs[NCOMP_MAX][NCOMP_MAX];
extern int Geqn_start[NCOMP_MAX];
extern int Nblock[NCOMP_MAX];
extern int Nseg_per_block[NCOMP_MAX][NBLOCK_MAX];
extern int SegType_per_block[NCOMP_MAX][NBLOCK_MAX];
extern int Grafted_Logical;
extern int Grafted[NCOMP_MAX];
extern int Graft_wall[NCOMP_MAX];
extern int GraftedWall_TF[NWALL_MAX_TYPE];
extern double *Poly_graft_dist;     /* distance associated with polymer grafting - */
extern double Rho_g[NCOMP_MAX];
extern double G_prefactor;
extern int Ntype_mer;
extern int Nmer[NCOMP_MAX];
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern int Npol_comp;
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
extern int Last_nz_cr;
extern int Nmer_t_total[NBLOCK_MAX];
extern int Type_mer_to_Pol[NBLOCK_MAX];
extern int Poly_to_Type[NCOMP_MAX][NBLOCK_MAX];
extern int Poly_to_Ntype[NCOMP_MAX];
extern int Nseg_tot;
extern int Nbond_max;
extern int Nseg_type[NCOMP_MAX];
extern int Icomp_to_polID[NCOMP_MAX];
extern int Grafted_SegID[NCOMP_MAX];
extern int Grafted_SegIDAll[NCOMP_MAX];
extern int Grafted_TypeID[NCOMP_MAX];
extern int **Nseg_type_pol;
extern char *Cr_file;
extern char *Cr_file2;
extern char cr_file_array[FILENAME_LENGTH];
extern char cr_file2_array[FILENAME_LENGTH];
extern char *Poly_file_name;
extern char poly_file_array[FILENAME_LENGTH]; 
extern char *WallPos_file_name;
extern char wallPos_file_array[FILENAME_LENGTH]; 
extern char *Outpath;
extern char Outpath_array[FILENAME_LENGTH]; 
extern char *Runpath;
extern char Runpath_array[FILENAME_LENGTH];
extern int Set_GUIDefaults_to_OLD_File;
extern int Read_OLDInput_File;
extern int Read_XMLInput_File;
extern char *InputOLD_File;
extern char InputOLDFile_array[FILENAME_LENGTH];
extern char EchoInputFile_array[FILENAME_LENGTH];
extern char *InputXML_File;
extern char InputXMLFile_array[FILENAME_LENGTH];
extern char *DensityFile;
extern char *DensityFile2;
extern char DensityFile_array[FILENAME_LENGTH];
extern char DensityFile2_array[FILENAME_LENGTH];
/*extern char Cr_file[FILENAME_LENGTH];
extern char Cr_file2[FILENAME_LENGTH];
extern char Cr_file3[FILENAME_LENGTH];
extern char Cr_file4[FILENAME_LENGTH];*/
extern double Cr_break[2];
extern int Ncr_files;
extern int SegAll_to_Poly[NMER_MAX];
extern int *Unk_to_Poly;
extern int *Unk_to_Seg;
extern int *Unk_to_Bond;
extern int ***Poly_to_Unk;
extern int **Poly_to_Unk_SegAll;
extern int Ngeqn_tot;
extern int Nbonds;
extern int **Nbond;
extern int ***Bonds;
extern int ***pol_sym_tmp;
extern int *Pol_Sym;
extern int *Pol_Sym_Seg;
extern int *BondAll_to_isegAll;
extern int *BondAll_to_ibond;
extern int Unk2Comp[NMER_MAX];
extern int Nmer_comp[NCOMP_MAX];
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
extern int **Bonds_SegAll;
extern int *Nbonds_SegAll;
extern double Gsum[NCOMP_MAX];	/* chain partition func for CMS_SCFT */
extern double *Total_area_graft;	/* Total surface area to use for grafted chains */
extern double *Gsum_graft;	/* prefactor term for grafted chains */
extern double *Gsum_graft_noVolume;	/* prefactor term for grafted chains */
extern double **GsumPrefac_XiDerivs; /* keep track of prefactors so we can implement Jacobians for tethered chains */
extern double ***GsumPrefac_GDerivs; /* keep track of prefactors so we can implement Jacobians for tethered chains */
extern int **Index_SurfNodes_Gsum; /* keep track of surface nodes we can implement Jacobians for tethered chains */
extern int ***Index_UnkGQ_Gsum; /* keep track of unknowns when computing Gsum for tethered chains */
extern int **Index_UnkB_Gsum; /* keep track of unknowns when computing Gsum for tethered chains */
extern int **Nodes_Surf_Gsum; /* counter for surface nodes used to compute Gsum for tethered chains */

/*some continuation related arrayes */
extern int  Cont_ID[NCONT_MAX][2];  /* Array of iwall/icomp ids for use in continuation.  */
extern int NID_Cont;

/*********************************************************************/
extern double Ads[NCOMP_MAX][2];
extern double Ads_ex[NCOMP_MAX][2];
extern double *Integration_profile; /* a place to put the integrand as a function of position */

/****************************************************************************/
/*#ifdef __cplusplus
}
#endif*/
