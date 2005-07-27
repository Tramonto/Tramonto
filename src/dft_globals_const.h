/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

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
#endif

#include "az_aztec.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"

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

/* basic constants */
#define PI    M_PI

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
#define NWALL_MAX_TYPE 50 

#define NBOND_MAX 4

#define NSTEPS_MAX 10

/* a constat flag to indicate that there is no bond between a pair of segments in 
   a polymer problems */
#define NO_BOND_PAIR -962.0

/* a constant flag for the share_global_vec routine */
#define UNINIT_VEC -200.0
/*
 * Stencil types refer to the integration schemes needed for different
 * non-local physics. 
 * NSTEN is the total number of stencil types.
 * The following are the current acceptable values for Stencil types.
 */

#define NSTEN        8

#define DELTA_FN       0
#define THETA_FN       1
#define U_ATTRACT      2
#define THETA_CHARGE   3
#define POLYMER_CR     4
#define POLYMER_GAUSS  5
#define THETA_FN_SIG   6
#define DELTA_FN_BOND  7

/*
 * An equation type list is given here to make identification of a particular
 * unknown number straightforward.  This became necessary with the introduction
 * of the WTC polymers where we now have 6 types of equations to fill.
 */
#define NEQ_TYPE       8 
#define NO_UNK        -888

#define DENSITY        0
#define RHOBAR_ROSEN   4
#define DIFFUSION      5
#define CAVITY_WTC     6
#define BOND_WTC       7

#define CMS_FIELD      1
#define CMS_G          2 
#define POISSON        3

/* Here are some constants needed to make the physics based ordering of the
   matrix an option in the code. */
#define MATRIX_FILL_NODAL 1   /* set to zero for physics based ordering */
#define NODAL_FLAG -999
#define BOX 0
#define LOCAL 1
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
#define cyl_periodic_3D                 6
#define cyl2D_sphere3D_pore             7
#define cyl3D_slit2D_pore               8
#define tapered_pore                    9

/* 
 * These constants identify the functional choices (Type_func).
 */
#define NONE       -1
#define ROSENFELD  0
#define ROSENFELD2 1
#define LDA        2
#define GHRM       3
#define GVDWM      4


/*
 * These constants identify attraction functional choices (Type_attr).
 */
#define NONE -1
#define LJ_WCA_CS 0

/* 
 * These constants identify the functional choices (Type_coul).
 */
#define NONE      -1
#define BARE       0
#define DELTAC     1 
#define POLARIZE   2
#define LIKE_LJ    3

/*
 *  These constants identify choice for polymers.
 */
#define NONE           -1
#define CMS_FR_JNT_CHN1     0
#define CMS_FR_JNT_CHN2 1
#define CMS_FR_JNT_CHN3 2
#define CMS_GAUSSIAN    3
#define CMS_SCFT        4
#define TC_FR_JNT_CHN   5

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

/*
 *  The following define ways to distribute charge in the system 
 */
#define POINT_CHARGE 0
#define SMEAR_CHARGE 1
#define BACKGROUND   2

/*
 * The following are choices for the initial guess (Iguess)
 */

#define CONST_RHO       -3 
#define CONST_RHO_L     -2 
#define CONST_RHO_V     -1 

#define EXP_RHO          0
#define EXP_RHO_L        1
#define EXP_RHO_V        2

#define STEP_PROFILE     3

#define CHOP_RHO         4
#define CHOP_RHO_L       5
#define CHOP_RHO_V       6
#define CHOP_RHO_STEP    7

#define LINEAR           8
/*
 * The following are the various fields for continuuation 
 */
#define CONT_MESH        0   /* mesh size */
#define CONT_TEMP        1   /* State Parameters */
#define CONT_RHO_0       2
#define CONT_RHO_ALL     3
#define CONT_LOG_RHO_0   4 
#define CONT_LOG_RHO_ALL 5 
#define CONT_SCALE_RHO   6

#define CONT_EPSW_0      7    /* Wall-Wall Energy Params */
#define CONT_EPSW_ALL    8
#define CONT_SCALE_EPSW  9

#define CONT_EPSWF00     10    /* Wall-Fluid Energy Params */
#define CONT_EPSWF_ALL_0 11 
#define CONT_SCALE_EPSWF 12

#define CONT_EPSFF_00    13   /* Fluid-Fluid Energy Params */
#define CONT_EPSFF_ALL   14   
#define CONT_SCALE_EPSFF 15

#define CONT_SCALE_CHG   16  /* Charged surface params */
#define CONT_SEMIPERM   17  /* Vext_membrane */
#define CONT_WALLPARAM  18  /* Vext_membrane */

#define CONT_CRFAC  19  /* continuous mixing of two cr files */


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
#define VEXT_1D          2  /* 1D potential for 1D systems */
#define VEXT_1D_XMIN     3  /* crude 1D-like treatment of funny geometries */
#define VEXT_1D_ORIENTATION     4  /* crude 1D-like treatment of funny geometries */
#define VEXT_3D_INTEGRATED      5  /* more proper 3D integration potential for funny geometries */
#define VEXT_ATOMIC        6  /* 3D potential for 3D problems */

 /* The following is a choice for the wall type Ipot_ww_n*/
#define NO_WW              0
#define ATOM_CENTERS_WW    1 

 /* The following is a choice for the BC type */
#define IN_WALL     -1
#define IN_BULK      0
#define PERIODIC     1
#define REFLECT      2
#define LAST_NODE    3

/* These choices are for types of boundary conditions for
    charged systems                                    */
#define CONST_POTENTIAL  1
#define CONST_CHARGE     2
#define ATOMIC_CHARGE     3

 /* The following is a choice for the load balance flag */
#define LB_LINEAR    0
#define LB_BOX       1
#define LB_WEIGHTS   2
#define LB_TIMINGS   3
#define LB_NON0      4
#define LB_MIXED     5

 /* The following is a choice for the Fast_fill_flag */
#define CHECK_NONE   0
#define CHECK_BC     1
#define CHECK_HW     2
#define CHECK_BOTH   3

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


/* The following are choices for output */
#define MINIMAL      0
#define DENSITIES    1 
#define NO_SCREEN    2 
#define VERBOSE      3 

/* The followint are choices for the output of density info */
#define SWITCH_RHO   0
#define SWITCH_RELP  1
#define SWITCH_ION   2
#define SWITCH_MU    3

#define SWITCH_SURFACE_SEP   0

/* The following is a flag for 1D boundary conditions in a 2D or
    3D domain --- only set up for steady state problems now */
#define FLAG_1DBC -999
#define FLAG_BULK -888
#define BULK_ZONE    2 

/* Polymer constants */
#define N_NZCR_MAX   200   /* maximum # of non-zero's in direct correlation fn */
#define NBLOCK_MAX   5
#define NMER_MAX     40


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

/*
 * Declaration of the Trilinos solver stuff follows.
 */
extern void * LinProbMgr_manager;

/* 
 *  Declaration of the Aztec_Struct follows. 
 *  Information for the linear solver is all held here.
 */

struct Aztec_Struct {
  int    options[AZ_OPTIONS_SIZE]; /* Array used to select solver options.  */
  double params[AZ_PARAMS_SIZE];   /* User selected solver paramters.       */
#ifdef DONE_WITH_THESE
  int    proc_config[AZ_PROC_SIZE];/* Processor information.                */
  int    *data_org;                /* Array to specify data layout          */
  double status[AZ_STATUS_SIZE];   /* Information returned from AZ_solve(). */
  int    *update;                  /* vector elements updated on this node. */
  int    *external;                /* vector elements needed by this node.  */
  int    *update_index;            /* ordering of update[] and external[]   */
  int    *extern_index;            /* locally on this processor.            */
  int    *bindx;                   /* Sparse matrix to be solved is stored  */
  double *val;                     /* in these MSR arrays.                  */
  int    N_update;                 /* # of unknowns updated on this node    */
  int    nonzeros;                 /* # of nonzeros in sparse matrix        */
#endif
};

/* 
 *  Declaration of the Loca_Struct follows.
 *  Information for the continuation library is held here.
 */

#define LOCA
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


/*extern
int    count_zero;
extern
int    count_nonzero;
*/

/* Basic Equation info */
extern int Phys2Nunk[NEQ_TYPE];  /* Number of unknowns of a particular equation type */
extern int Phys2Unk_first[NEQ_TYPE]; /* starting unknown number for a given equation type */
extern int Phys2Unk_last[NEQ_TYPE]; /* ending unknown number for a given equation type */
extern int Unk2Phys[3*NCOMP_MAX+NMER_MAX+NMER_MAX*NMER_MAX+13]; /* array that gives equation type
                                                                         given an unknown index */

/* Mesh info */


/*************** Global Mesh ********************************/

extern
int     Ndim;            /* # of spatial dimensions of the problem      */
extern
int     Nnodes;          /* # of nodes in the mesh                */
extern
int     Nunknowns;       /* # of unknowns in the problem          */
extern
int     Nelements;       /* # of elements in the mesh             */
extern
int     Elements_plane;  /* # of elements in the x1-x2 plane      */
extern
int     Nodes_plane;     /* # of nodes in the x1-x2 plane          */
extern
int Nodes_x[NDIM_MAX];   /* Array[Ndim]: # nodes in each dimension */
extern
int Elements_x[NDIM_MAX];/* Array[Ndim]: # elements in each dim    */
extern
int     Max_sten_length[3];  /* The number of nodes in the longest stencil */
extern
int     Sten_length_hs[3];  /* The number of nodes in the h.s. stencil */

/*************** Reference Variables ***********************/
extern
double Length_ref;
extern
double Density_ref;
extern
double Dielec_ref;
extern
double Potential_ref;
extern
double VEXT_MAX;

/*************** Extended Local Mesh ************************/
extern
int     Nnodes_box;          /* # of nodes in the extended local mesh   */
extern
int     Nunknowns_box;       /* # of unknowns in the extended local mesh */
extern
int     Nelements_box;       /* # of elements in the extended local mesh */
extern
int     Elements_plane_box;  /* # of el in x1-x2 plane  of extended local mesh  */
extern
int     Nodes_plane_box;     /* # of nodes in x1-x2 plane of extended local mesh*/
extern
int Nodes_x_box[NDIM_MAX];   /* Array[Ndim]: # nodes in idim on extended local mesh */
extern
int Elements_x_box[NDIM_MAX];/* Array[Ndim]: # elements in idim on extended local mesh */
extern
int     Min_IJK_box[3]; /* The minimum IJK of the extended internals on this procesor */
extern
int     Max_IJK_box[3];  /* The maximum IJK of the extended internals on this procesor */
extern
int     Pflag[3];   /* PERIODIC flag array[Ndim]: TRUE Nodes_x_box=Nodes_x */

/************** Local Mesh **********************************/
extern
int     Min_IJK[3];      /* The minimum IJK of the internals on this procesor */
extern
int     Max_IJK[3];      /* The maximum IJK of the internals on this procesor */

/************** Communications arrays for Gather vectors ******/
extern
int *Comm_node_proc;   /* array on proc 0 of nodes per processors */
extern
int *Comm_unk_proc;   /* array on proc 0 of unknowns per processors */
extern
int *Comm_offset_node; /* array on proc 0 of offsets of nodes  MPI_Gatherv*/
extern
int *Comm_offset_unk;  /* array on proc 0 of offsets of unknowns MPI_Gatherv*/

/************** Other **************************************/
extern
double  Vol_el;          /* Volume of one element of our regular mesh        */
extern
double  Area_surf_el[3]; /*Area of surface element with normal in idim direction*/
extern
double Vol_in_surfs[NCOMP_MAX];  /* volume in all of the surfaces for each list */
extern  
int     Nlists_HW;       /* Number of lists needed if hard walls (for mixtures)*/
extern
int     Nel_wall_count;  /* Number of elements in the 0.5sigma of wall:rhobars*/
extern
int    **Nelems_per_wall;   /* Number of elements in a given [iwall] wall        */
extern
int     Nnodes_per_proc;  /* Number of nodes owned by this processor         */
extern
int     Nunk_int_and_ext; /* Number of unknownsneeded on this processor      */
extern
int     *B2L_unknowns;     /* Box to Local array for all unknowns */
extern
int     *B2G_node;         /* Box to global array for all box nodes */
extern
int     *B2G_unk;          /* Box to global array for all box unknowns */
extern
int     *L2B_node;         /* Local to box array for all local nodes */
extern
int     *B2L_node;         /* Box to local array for all local nodes */
extern
int     *L2G_node;         /* Local to global coordinates */
extern
int     Nunk_per_node;   /* Number of unknowns per node (usually Ncomp       */
extern
int     Nrho_bar;        /* Number of rhobar equations per node */
extern
int     Nrho_bar_s;      /* Number of scalar rhobar equations per node */
extern
int     Npoisson;        /* Number of rhobar equations per node */
extern
int     Ndiffusion;        /* Number of rhobar equations per node */
extern
int     Ndensity_unk;         /* Number of unknowns for the Euler-Lagrange equation */
extern
int     Ntype_unk;       /* Number of equations defining segment types for polymer TC cases */
extern
int     Nrho_bar_cavity; /* Number of nonlocal densities for the cavity function - WTC polymers */
extern
int     Nrho_bar_bond; /* Number of nonlocal densities for the bond functionals - WTC polymers */
extern
int     Nnodes_per_el_V;  /* Number of nodes per volume element              */
extern
int     Nnodes_per_el_S;  /* Number of nodes per surface element            */
extern
int     Plane_new_nodes; /* Indicates in which plane (xy,yz,xz)nodes are added*/
extern
int     Pos_new_nodes;   /* Indicates where nodes are added                  */
extern
double  Size_x[NDIM_MAX];    /*Array of the size of the domain in each dim. */
extern
double  Esize_x[NDIM_MAX];   /*Array of the size of an element in each dim. */
extern
int     Lmesh_refine;       /*Switch for auto mesh refinement               */
extern
int     Type_bc[NDIM_MAX][2];/*Array of boundary conditions in each dim.    */
extern
int     Non_unique_G2B[4];   /* Flag indicating which box dimensions       
                                wrap around completely, resulting in a
                                non-unique G2B mapping */
extern
int   **Nodes_2_boundary_wall;  /*Array[Nlists_HW][Nnodes] -1 if nod b.n. else b.n.  */
extern
int   **Wall_elems;     /*Array[Nlists_HW][Nelements] TRUE for wall elements */
extern
int   ****Touch_domain_boundary; /*Array[Nwall][Nlists_HW][[Ndim][2] =0 if surface hits left
                                  boundary = 1 if hits right domain boundary else -1 */
extern
int   **Xtest_reflect_TF; /* Array[Nwall][Ndim] for reflections/wall-wall boundaries */
/* these are used to set up the Nodes_2_boundary_wall array then discarded*/
extern
int   Nnodes_wall_box; /* Count number of nodes in box that touch a wall */
extern
int  *Nodes_wall_box; /* Array to store which nodes touch a wall */
extern
int  *Nwall_touch_node; /* Array to store number of walls touching a given node */
extern
int  **Wall_touch_node; /*Array to store which walls touch a given node */
extern
int  **List_wall_node; /*Array to store which walls touch a given node */
extern
int  *Index_wall_nodes; /* ArraY to store indexing in these mesh arrays */

extern
int First_time; /* for MSR preprocessing */

extern
int     Nzone;          /* Number of diff. quadrature zones on the mesh      */
extern
int    *Nodes_to_zone;   /* Array[Nnodes] of quadrature zones */
extern
double  Rmax_zone[5];    /* Array distances from surfaces in quadrature zones */
extern
int     Mesh_coarsening;  /* Flag indicating whether mesh coarsening is on */
extern
int    *Mesh_coarsen_flag;/* Flag (Nnodes) telling how much coarsening for
                             a given node, or negative values telling which*/
extern
int     Coarser_jac;     /* Flag to switch on coarser jacobian than residual */
extern
double  Jac_grid;     /* Flag to switch on coarser jacobian than residual */
extern
int     Lcut_jac;  /* Logical to indicate if Jacobian stencils will be cut off */
extern
double  Jac_threshold; /* Threshold level for Jacobian stencils ... max/Jac_threshold */
extern  
int    *Fast_fill_TF;    /* Flag for fast Jacobian fill */

/* Continuation info */
extern
int     Nodes_old;         /* # of nodes in the mesh of the previous run */
extern
int Nodes_x_old[NDIM_MAX];/* Array[Ndim]: # nodes in idim of previous run */
extern
double  *X_old;           /* Array of unknowns from previous run */
extern
double  *X2_old;           /* Array of unknowns from previous run */

extern
int     Print_force_type;  /* flag for printing of force */
extern
int     Print_rho_type;  /* flag for file type for printing of densities */
extern
int     Print_rho_switch; /* flag for printing densities -- format */
extern
int     Lprint_gofr; /* flag for printing radial distribution functions */
extern
int     Lprint_pmf; /* flag for printing radial distribution functions */
extern
int     Print_mesh_switch; /* flag for printing densities -- format */
extern
int     Lper_area;  /*logical for per unit are outputs of params */ 
extern
int     Lcount_reflect;  /*logical for per unit are outputs of params */ 
extern
int     Nruns;           /* Number of runs to perform (varying the mesh)     */
extern
double  Del_1[NWALL_MAX_TYPE];    /*Stepping parameter for field #1  */
extern
double  Guess_range[2];   /* surf sep to switch between using Rho_b and X_old */
extern
double  Rho_max;         /* max rho when using an old solution for mesh contin*/
extern
int     Imain_loop;    /* Couter on the number of times through the program  */
extern
double  Scale_fac; /* for continuation in arrays !*/


/* Surface Physics info */
extern
int     Nwall;           /* Number of surfaces in the calculation            */
extern 
int     Nwall_type;       /* Number of surface types in the problem */
extern
int     WallType[NWALL_MAX]; /* array containing type number for each surface */
extern
int     Nlink;           /* Number of macro-surfaces */
extern
int     Link[NWALL_MAX]; /* Index iwall to a particular macro-surface */
extern
int     *Nwall_this_link; /* Number of walls linked together in a particular macro-surface */
extern
int     **Link_list; /* List of walls linked in a particular macro-surface */
extern
int     Surface_type[NWALL_MAX_TYPE];    /* Type of surfaces of interest                     */
extern
int     Orientation[NWALL_MAX_TYPE];  /* Orientation of planar/bumpy infinite walls*/
extern
double  WallParam[NWALL_MAX_TYPE];/* Array[Nwall] of a characteristic wall parameter*/
extern
double  WallParam_2[NWALL_MAX_TYPE];/* Array[Nwall] of a characteristic wall parameter*/
extern
double  WallParam_3[NWALL_MAX_TYPE];/* Array[Nwall] of a characteristic wall parameter*/
extern
double  WallPos[NDIM_MAX][NWALL_MAX]; /* Array of the centers of the surfaces*/

/* Fluid Physics info */
extern
int     Ncomp;           /* Number of components in the current problem      */
extern
double  Temp;            /* Reduced Temperature for LJ interactions          */
extern
double  Temp_elec;       /* Reduced Temperature for Poisson's equation       */
extern
double  charge_fluid;    /* the charge in the fluid ... post processing */
extern  
int     Type_dielec;      /* choose how to handle dielectric constants in system */
extern
double  Dielec_bulk;     /* the dielectric constant in the bulk fluid */
extern
double  Dielec_pore;     /* the dielectric constant in the "pore" fluid */
extern
double  Dielec_X;        /* distance from a surface defines the "pore" fluid */
extern
double *Dielec_wall;     /* the dielectric constant as function of wall type */
extern
double *Dielec;     /* dielectric constant as function of ielement_box */
extern
int     Lpolarize[NCOMP_MAX]; /* logical for if a species is polarizeable */
extern
double  Pol[NCOMP_MAX];  /* bulk polarizeability for each species */
extern
double  Energy;   /* Surface free energy to return to Towhee */


extern
double  Betap;           /* Pressure in units of kT sigma_ff[1]^3           */
extern
double Betap_LBB;       /* Pressure calculated for LBB of domain */
extern
double Betap_RTF;       /* Pressure calculated for RTF of domain */
extern
double  Betap_id;       /* Ideal gas Presseure in units of kT sigma_ff[1]^3   */
extern
double  Betap_att;      /* Attractive Presseure in units of kT sigma_ff[1]^3   */
extern
double  P_over_po;
extern
double  Hs_diam[NCOMP_MAX]; /* Array of effective hard sphere diameters      */
extern
int     L_isotherm; /* Logical for isotherm vs. force per distance data */
extern
double  Rho_b[NCOMP_MAX];   /* Array[Ncomp] of component bulk densities      */
extern
double  Rhobar_b[10]; /* Array[Nrho_bar] of bulk rhobars      */
extern
double  Rhobar_b_LBB[10]; /* Array[Nrho_bar] of bulk rhobars LBB  */
extern
double  Rhobar_b_RTF[10]; /* Array[Nrho_bar] of bulk rhobars RTF  */
extern
double  Rho_seg_b[NMER_MAX]; /* array of bulk segment densities */
extern
double  Rho_seg_LBB[NMER_MAX];
extern
double  Rho_seg_RTF[NMER_MAX];
extern
double  *Rhobar3_old;   /* Array[Nnodes_box] of old values of rhobar 3*/
extern
double Xi_cav_b[4]; /* Array of bulk rhobars for cavity functions of WTC polymer functionals */
extern
double Xi_cav_LBB[4]; /* Array of bulk rhobars for cavity functions of WTC polymer functionals */
extern
double Xi_cav_RTF[4]; /* Array of bulk rhobars for cavity functions of WTC polymer functionals */
extern
double BondWTC_b[NMER_MAX*NMER_MAX]; /*Array of bulk rhobars for bonds in WTC functionals*/
extern
double BondWTC_LBB[NMER_MAX*NMER_MAX]; /*Array of bulk rhobars for bonds in WTC functionals*/
extern
double BondWTC_RTF[NMER_MAX*NMER_MAX]; /*Array of bulk rhobars for bonds in WTC functionals*/
extern
double  Rho_coex[2];   /* Liquid and Vapor Coexisting Densities         */
extern
double  Betamu_hs_ex[NCOMP_MAX];/* Array of excess hardsphere chemical potentials*/
extern
double  Betamu[NCOMP_MAX];   /* Array[Ncomp] of chemical potentials*/
extern
double  Betamu_id[NCOMP_MAX];   /* Array[Ncomp] of ideal gas chemical potentials*/
extern
double Betamu_wtc[NMER_MAX];
extern
double  Betamu_ex_bondTC[NCOMP_MAX][NMER_MAX*NMER_MAX];/* Array of excess segment chemical potentials - WTC poolymer*/
extern
double  Betamu_seg[NMER_MAX];/* Array of excess segment chemical potentials - WTC poolymer*/
extern
int     Ipot_ff_n;    /* Potential Type for neutral part of f-f interactions */
extern
int     Ipot_wf_n[NWALL_MAX_TYPE];    /* Potential Type for neutral part of w-f interactions */
extern
int     Ipot_ww_n[NWALL_MAX_TYPE][NWALL_MAX_TYPE];    /* Potential Type for neutral part of w-f interactions */
extern
int     Ipot_ff_c;    /* Potential Type for charged part of f-f interactions */
extern
int     Ipot_wf_c;    /* Potential Type for charged part of w-f interactions */
extern 
int     Lhard_surf;   /* Logical indicating if the surfaces have hard cores */
extern  
int     Iliq_vap;     /* Type of liquid vapor profile */
extern
int     Iguess1;        /* Type of initial guess */
extern
int     Nsteps;         /* Number of steps for a step profile initial guess */
extern
int     Orientation_step[NSTEPS_MAX]; /* orientation of the step profile */
extern
double  Xstart_step[NSTEPS_MAX];  /* start position array for the step profile */
extern
double  Xend_step[NSTEPS_MAX];  /* end position array for the step profile */
extern
double  Rho_step[NCOMP_MAX][NSTEPS_MAX];  /* density array for a step profile */
extern
int     Lbinodal;        /* Logical TF for binodal calculation */
extern
double  Thickness;    /* Thickness parameter for doing wetting studies */

extern
double  Alpha;        /* Yukawa decay parameter                              */
extern
int     Mix_type;  /* Choice of Mixing Rules */
extern
double  Mass[NCOMP_MAX];           /* Array of the mass of each specie*/
extern
double  Sigma_ff[NCOMP_MAX][NCOMP_MAX];/* Array of f-f interaction diameters */
extern
double  Bond_ff[NCOMP_MAX][NCOMP_MAX];/* Array of f-f bond lengths for polymers */
extern
double  Eps_ff[NCOMP_MAX][NCOMP_MAX];  /* Array of f-f interaction energies  */
extern
double  Cut_ff[NCOMP_MAX][NCOMP_MAX];  /* Array of f-f cutoff distances      */
extern
double  Charge_f[NCOMP_MAX];           /* Array of the valence of each specie*/
extern
double  Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];/* Array of w-f interaction diameters */
extern
double  Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];  /* Array of w-f interaction energies  */
extern
double  Cut_wf[NCOMP_MAX][NWALL_MAX_TYPE];  /* Array of w-f cutoff distances      */
extern
double  Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];/* Array of w-w interaction diameters */
extern
double  Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];  /* Array of w-w interaction energies  */
extern
double  Cut_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];  /* Array of w-w cutoff distances      */

extern
int     **Lsemiperm; /*Array of logicals for semi-permeable surfaces */
extern
double  **Vext_membrane; /*Array potentials for semi-perm surfaces */
extern
double  **Vext_set;      /*Array of maximum set points for ext potentials */
extern
double **Vext;        /* External field array [Nnodes][Ncomp]           */
extern
double **Vext_coul;        /* Coulomb External field array [Nnodes]           */
extern
double *Vext_old;        /* for post processing: ext field array           */
extern
double ***Vext_dash;  /* derivative of external field [Nnodes][Ncomp][Nwall]*/
extern
double **Uww;        /* Wall-Wall interactions [Nwall-1][Nwall-1]           */
extern
double **Uww_link;        /* Wall-Wall interactions [Nlink-1][Nlink-1]           */
extern
double **X_wall;  /* Distance from inode to iwall [Nnodes][Nwall]    */
extern
double **X_wall2;  /* Distance from inode to iwall [Nnodes][Nwall]    */
extern
int    **Zero_density_TF; /* array [Nnodes][icomp] indicates where VEXT_MAX */
extern
double  *Betamu_att;   /* sum over jcomp of Van der waals constant a(icomp,jcomp)*/
extern
double  **Avdw;    /*  Van der waals constant a(icomp,jcomp)*/

extern
double  Sigma_w[NWALL_MAX_TYPE];  /* Array[Nwall] of w-w interaction diameters    */
extern
double  Eps_w[NWALL_MAX_TYPE];    /* Array[Nwall] of w-w interaction energies     */
extern
double  Rho_w[NWALL_MAX_TYPE];    /* Array[Nwall] of w-w interaction energies     */
extern
double  Elec_param_w[NWALL_MAX]; /* Array: surf charge(potential) per area  */
extern
int     Type_bc_elec[NWALL_MAX_TYPE];/* Array of surface b.c.'s for charged systems    */
extern 
int     Nlocal_charge; /*Number of localized surface charges in the system */
extern
double  **Charge_x;    /*Position of a given local charge [Nlocal_charge][Ndim]*/
extern
double  *Charge_Diam;  /*Diameter of a given local charge [Nlocal_charge]*/
extern
double  *Charge;  /*Value of the local charge [Nlocal_charge]*/
extern
int     Charge_type_atoms; /* Type of charge distribution on the atoms */
extern
int     Charge_type_local; /* Type of charge distribution on the added local charge */
extern
double  *Deltac_b;   /* Array [icomp] of electrostatic correlations in bulk*/
extern
double  **Charge_w_sum_els; /*Array[Nnodes_b][Ndim] of surface charge per area*/
extern
double  *Charge_vol_els; /*Array[Nelemts_box] of volume charge per element */
extern
int     Vol_charge_flag; /* Flag for volumetric charges */
extern
int     Surf_charge_flag; /* Flag for volumetric charges */
extern
int     **Nelems_S; /* Array[Nlists_HW][Nnode_per_proc] of # surf elems the b.node touches*/
extern
int     ***Surf_normal; /*Array[Nlists][Nnodes_per_proc][Nelems_S] of unit normal vectors */
extern
int     ***Surf_elem_to_wall; /*Array of wall to which a given surface element belongs...
                                 a given node may belong to more than one wall !! */
extern
int     **Surf_elem_type; /*Array[Nnodes_per_proc][Nelems_S] of surf elem type */
extern
double  ***S_area; /*Array[Nlists][Nwall][Ndim] of total surf. area in idim on iwall*/ 
extern
double  **S_area_tot; /*Array[Nlists][Nwall] total surf. area on iwall*/ 
extern
int     **Nwall_owners; /*Array[Nilists][el_box] of number of walls 
                         (including images) that own a given element */
extern
int     ***Wall_owners; /*Array[ilist][iel_box][Nwall_owners] that stores
                          all of the wall owners of a given element */ 

/* Steady State Solutions Info */
extern
int    Lsteady_state;          /*True-False Steady State or Equilibrium Run*/
extern
int    Linear_transport;       /*True-False Steady State or Equilibrium Run*/
extern
double Velocity;               /*Constant Convective Velocity over Diffusion coefficient*/
extern
int    Grad_dim;               /*direction where gradient is implemented*/
extern
int    L1D_bc;               /*logical for 1D boundary condition in Grad_dim direction */
extern
double X_1D_bc;               /*distance where 1D boundary is applied */
extern
double X_const_mu;             /*length where constant mu applies @ edges of domain*/
extern
double Rho_b_LBB[NCOMP_MAX];   /*Rho_b boundary condition left-bottom-back*/
extern
double Rho_b_RTF[NCOMP_MAX];   /*Rho_b boundary condition right-top-front*/
extern
double Elec_pot_LBB;           /*Electric potential boundary condition LBB*/
extern
double Elec_pot_RTF;           /*Electric potential boundary condition RTF */
extern
double Betamu_LBB[NCOMP_MAX];  /*Chemical Potential Boundary Condition LBB */
extern
double Betamu_RTF[NCOMP_MAX];  /*Chemical Potential Boudary Condition RTF*/
extern
double D_coef[NCOMP_MAX];  /*Diffusion Coefficients for ion species */
extern
double *Pore_rad_L_IC;    /* array of left  Radii of ion chan pore segments (1D) */
extern
double *Pore_rad_R_IC;    /* array of right Radii of ion chan pore segments(1D) */
extern
double *Lseg_IC;          /* array of length of ion chan pore segments (1D) */
extern
int    Nseg_IC;           /* number of pore segments in a given ion channel */
extern
double    *Area_IC;      /* 1D ion channel area per node (box units)*/
extern
int   Geom_flag;    /* geometry flag for ion chan. see OPTION_ definitions*/

/* OUTPUT INTEGRAL PARAMETERS */
extern
int    ***Nel_hit;      /* number of elements hit by a given node in a given list */
extern
int    ***Nel_hit2;     /* same as prev. for a bulk fluid */
extern
int    List[2];       /* which list numbers we care about for integrals*/
extern
int    Imax;          /* how many lists are relevent to the case at hand 1 or 2 */

/* SOME CONSTANTS */
extern
double Inv_4pi;               /* Precalculated value of 1/(4*pi)                    */
extern
double Inv_rad[NCOMP_MAX];    /* Precalculated inverse of component radius          */
extern
double Inv_4pir[NCOMP_MAX];   /* Precalculated inverse of component's 4 pi Radius   */
extern
double Inv_4pirsq[NCOMP_MAX]; /* Precalculated inverse of component's 4 pi Radius^2 */

extern
int     Type_func;    /* Type of functional for the calculation              */
extern
int     Type_attr;    /* Type for handling attractions                       */
extern
int     Type_coul;    /* Type for handling coulomb interactions              */
extern
int     Type_poly;    /* Type for handling polymers                          */
extern
int     Type_poly_TC;  /* Type for Wertheim-Tripathi-Chapman polymers */

/* Startup Info */
extern
int     Restart;     /* Logical that switches between new prof & restart file*/
extern
int     Restart_field[NEQ_TYPE];
extern
int     Iwrite;       /* Do we want a complete or modified set of output data*/

/* Parallel Info, Aztec info */
extern
int     Num_Proc; /* The total number of processors used in this calculation */
extern
int     Proc;     /* The unique  processor number (from 0 to Num_Proc-1)     */
extern
struct  Aztec_Struct Aztec; /* Structure to hold all the Aztec info          */
extern
int     Load_Bal_Flag; /* Flag specifying type of laod balancing to do       */


/* Nonlinear Solver info */
extern
int Max_Newton_iter;    /* Maximum # of Newton iterations (10 - 30)          */
extern
double Newton_abs_tol,Newton_rel_tol; /* Convergence tolerances (update_soln)*/
extern
double Min_update_frac; /* Minimum fraction to update solution to slow down
                           Newton's method */

/* Linear Solver info */
extern
int Az_solver;
extern
int Az_kspace;
extern
int Az_scaling;
extern
int Az_preconditioner;
extern
double Az_ilut_fill_param;
extern
double Az_tolerance;
extern
int    Max_gmres_iter;


extern
struct  Loca_Struct Loca; /* Information for continuation library */
/*
 * The global variable Stencil is a 3D array of Stencil_Struct
 * (for hard spheres of size [Nsten][Ncomp][Ncomp] ). This variable is
 * defined in dft_stencil.h, and the extern statement is found here.
 */

extern 
struct Stencil_Struct ***Stencil;
extern                    
int MPsten_Npts_R[NZONE_MAX];  /* # radial gauss pts. in MIDPOINT rule */
extern                    
int MPsten_Npts_arc[NZONE_MAX]; /* # theta gauss pts. in MIDPOINT rule */
extern                    
int MPsten_Npts_phi[NZONE_MAX]; /* # phi gauss pts. in MIDPOINT rule */
extern 
int Sten_Type[NSTEN]; /* on/off Flag for stencil types  */
extern
int Sten_Choice_S[NSTEN][NZONE_MAX]; /* Quadrature type fore each stencil   */
extern
int Sten_Choice_R[NSTEN][NZONE_MAX]; /* # Radial Gauss points for THETA_FNCs */

/* Polymer variables */
extern
double Deltar_cr,Gauss_a,Gauss_k,***Rism_cr;
extern
double Crfac;
extern
double Cr_rad[NCOMP_MAX][NCOMP_MAX];
extern
double Cr_rad_hs[NCOMP_MAX][NCOMP_MAX];
extern 
double Bupdate_fact;
extern
int Bupdate_iters,Geqn_start[NCOMP_MAX];
extern
int Nblock[NCOMP_MAX],Ntype_mer,Nmer[NCOMP_MAX],Type_mer[NCOMP_MAX][NMER_MAX];
extern
int Npol_comp,Nmer_t[NCOMP_MAX][NBLOCK_MAX],Last_nz_cr;
extern int Nmer_t_total[NBLOCK_MAX];
extern int Nseg_tot;
extern int Nseg_type[NCOMP_MAX];
extern char Cr_file[40],Cr_file2[40],Cr_file3[40],Cr_file4[40];
extern double Cr_break[2];
extern int Ncr_files;
extern int *Unk_to_Poly, *Unk_to_Seg, *Unk_to_Bond;
extern int ***Poly_to_Unk,**Poly_to_Unk_SegAll;
extern int Ngeqn_tot, Nbonds, **Nbond,***Bonds;
extern int *Pol_Sym;
extern int Unk2Comp[NMER_MAX],SegChain2SegAll[NCOMP_MAX][NMER_MAX],**Bonds_SegAll,*Nbonds_SegAll;


/*********************************************************************/
extern
double Ads[NCOMP_MAX][2];
extern
double Ads_ex[NCOMP_MAX][2];

/****************************************************************************/

/*** prototypes for function calls to another file  ***/
extern double second(void);

extern void read_input_file(char *,char *);
extern void thermodynamics(char *, int);
extern double dmu_drho_hs(double *);
extern double dmu_drho_att(double *);
extern void calc_stencils(void);
extern void set_gauss_quad(int , double *, double *);
extern void set_up_mesh(char *,char *);
extern void load_balance(int, double *, int *, int **);
extern void boundary_setup(char *);
extern void boundary_free(void);
extern void control_mesh(FILE *, char *,int, int *);
extern void free_mesh_arrays(void);
extern void setup_external_field_n(int **, int ***);
extern void  setup_vext_coulomb_vol();
extern void read_external_field_n(char *);
extern void read_zero_density_TF(char *);
extern void set_initial_guess(int,double **);
extern int  solve_problem(double **, double **);
extern int  newton_solver(double **x, void *con_ptr);

/*EXTERNAL FIELD AND WALL_WALL INTERACTION ROUTINES */
extern double integrate_potential(int, double, double, double, int, int, 
                      double *, double *, double *, double *,
                      double *, double *);

extern double uLJ12_6_cut(double,double, double, double);
extern double uderiv_LJ12_6(double,double,double,double,double);
extern double uLJ_wp(double,int,int);
extern double Vext_1D(double,int,int);
extern double Vext_1D_dash(double,int,int);
extern double uCOULOMB(double,double,double);
extern void   setup_wall_wall_potentials(int **, int ***);
extern void find_images(int, double, int *, double **,
                 double *, double *);
extern void find_images_1D(int, double, int *, double **,
                 double *, double *);
extern void find_images2(int, double, int *, double **,
                  double *, int,int);
extern void find_images_coulomb(int,int *, double **, double *);




/* RESIDUAL AND MATRIX FILL ROUTINES */
extern int loc_find(int,int,int);
extern void fill_resid_and_matrix_control(double **, int, int);
extern void fill_resid_and_matrix(double **, int, int, int);
extern void fill_resid_and_matrix_rb(double *, double *, 
                                  int **, double *, int, int);
extern void fill_resid_and_matrix_P(double **, int, int, int);
extern int  get_integration_pts(int, int, double ***, double **);
extern void pre_calc_rho_bar(struct RB_Struct *, double *, int,int,
                             double ***, int ***,double *);
extern void pre_calc_coarse_rho_bar(struct RB_Struct *);
extern void  pre_calc_dphi_drb(struct RB_Struct *, struct RB_Struct *,
                        struct RB_Struct *, struct RB_Struct *,
                        struct RB_Struct *);
extern void  pre_calc_dphi_drb1(struct RB_Struct *, struct RB_Struct *,
                        struct RB_Struct *, struct RB_Struct *,
                        struct RB_Struct *);
extern void  pre_calc_dphi_drb2(struct RB_Struct *, struct RB_Struct *,
                        struct RB_Struct *, struct RB_Struct *,
                        struct RB_Struct *);
extern void  pre_calc_dphi_drb_rb(struct RB_Struct *, double *,
                        struct RB_Struct *, struct RB_Struct *,
                        struct RB_Struct *, struct RB_Struct *,
                        double *);

extern void  pre_calc_dphi_drb_rb1(struct RB_Struct *, double **,
                        struct RB_Struct *, struct RB_Struct *,
                        struct RB_Struct *);

extern void  pre_calc_dphi_drb_rb2(struct RB_Struct *, double **,
                        struct RB_Struct *, struct RB_Struct *,
                        struct RB_Struct *);

extern double load_nonlocal_hs_rosen_rb(int, int, int,int,int, int,
                       int *,double **,
                       struct RB_Struct *,struct RB_Struct *,
                       struct RB_Struct *,struct RB_Struct *, int);

extern double load_rho_bar_s(int,double **,int, int,int,int,int *,int);

extern double load_rho_bar_v(double **,int,int,int,int,int *,int);

extern double load_mean_field(int, int, int, int, int, int *, double **, int);

extern double load_poissons_eqn(int, int, int, int *, double **);
extern double load_polarize_poissons_eqn(int, int, int, int *, double **);
extern double load_nonlinear_transport_eqn(int,int,int,int *, double **);
extern double load_linear_transport_eqn(int,int,int,int *, double **);
extern double load_poisson_bc(int,int,int);
extern double load_cavity_wtc(int,int,int,int,int *,double **);
extern double load_bond_wtc(int,int,int,int,int *,double **);
extern double load_polyTC_diagEL(int,int,int,int,int,int *,double **,int);
extern double load_polyTC_bondEL(int,int,int,int,int,int *,double **,int);
extern double load_polyTC_cavityEL(int,int,int,int,int,int *,double **,int);
extern void setup_polymer_cr(void);

extern double constant_boundary(int, int);
extern int find_jzone(int);

/*  MESH TRANSLATION ROUTINES */
extern int  map_0th_plane(int, int);
extern int  offset_to_node(int *,int *, int *);
extern int  offset_to_node_box(int *,int *, int *);
extern void node_to_ijk(int, int *);
extern int  ijk_to_node(int *);
extern void node_box_to_ijk_box(int, int *);
extern int  ijk_box_to_node_box(int *);
extern void ijk_to_ijk_box(int *,int*);
extern void ijk_box_to_ijk(int *,int*);
extern int  node_box_to_node(int);
extern int  unk_box_to_unk(int);
extern int  el_box_to_el(int);
extern int  el_to_el_box(int);
extern int  node_to_node_box(int);
extern int  node_to_node_box_no_bound(int);
extern int  unk_to_unk_box(int);
extern int  round_to_int(double);
extern void node_to_position(int,double *);
extern int  element_to_node(int);
extern int  element_box_to_node_box(int);
extern void element_to_nodes(int,int *);
extern void boundary_condition (int *);
extern void setup_surface(FILE *,int *, int **, int **, int ***, int *,int*,int ***);
extern int  node_to_elem(int, int, int *);
extern int  node_to_elem_return_dim(int, int, int *,int *, int *,int *);
extern int  node_box_to_elem_box_reflect(int, int, int *);
extern int  node_to_elem_v2(int, int);
extern void rescale_vext_max(void);
extern double HW_boundary_weight(int,int, double *, int, int *);
extern void find_match_dim(int,int, int **,int *,int ***);

/* SETUP THERMO ROUTINES */
extern void pot_parameters(char *);
extern double uLJatt_n(double,int,int);
extern double uLJatt_n_int(double,int, int);
extern double uLJatt_n_noshift(double,int, int);
extern double uCOULOMB_att(double,int, int);
extern double uCOULOMB_att_int(double,int, int);
extern double deltaC_MSA(double,int,int);
extern double deltaC_MSA_int(double,int,int);
extern double y_cav(double,double,double,double);
extern double dy_dxi2_cav(double,double,double,double);
extern double dy_dxi3_cav(double,double,double,double);



/* COMMUNICATIONS ROUTINES */
extern double gsum_double(double);
extern double gmax_double(double);
extern double gmin_double(double);
extern int    gsum_int(int);
extern int    gmax_int(int);
extern int    gmin_int(int);
extern void   gsum_int_vec(int *, int *, int);

/* POSTPROCESSING ROUTINES*/
extern void post_process(double **,char *,int *,double *,int, int);
extern double calc_adsorption(FILE *, double **,double,double);
extern void calc_surface_charge(FILE *, double **,double,double);
extern double calc_free_energy(FILE *, double **,double,double, int);
extern void calc_flux(FILE *,char *,double *);
extern void calc_force(FILE *, double **,double);
extern double calc_free_energy_polymer(FILE *,double **,double,double);

/* PRINTING - OUTPUT ROUTINES */
extern void collect_x_old(double **);
extern void collect_vext_old(void);
extern void print_profile_box(double **,char *);
extern void print_profile(char *);
extern void print_gofr(char *);
extern void print_vext(double **, char *);
extern void print_zeroTF(int **, char *);
extern void print_charge_vol(double *, char *);
extern void print_charge_surf(double **, char *);
extern void print_freen_profile_1D(double *, char *);
extern void print_Nodes_to_zone(int *, char *);
extern void print_time_histogram(int *,int *);
extern void setup_integrals(void);
/****************************************************************************/

/* 
   read all prototypes from the dft_*.c files to verify 
   that all routines are called with the correct signature
*/
#ifdef AGSGETRIDOFTHIS
  #include "include/dft_fill_shared.h"
  #include "include/dft_out_energy.h"
#ifdef PARALLEL
  #include "include/dft_comm.h"
#endif
  #include "include/dft_guess.h"
  #include "include/dft_out_flux.h"
  #include "include/dft_continuation.h"
  #include "include/dft_input.h"
  #include "include/dft_out_force.h"
  #include "include/dft_fill_main.h"
#ifdef RCB_LOAD_BALANCE
  #include "include/dft_ldbal.h"
#endif
  #include "include/dft_out_main.h"
  #include "include/dft_fill_mf.h"
  #include "include/dft_main.h"
  #include "include/dft_out_profiles.h"
  #include "include/dft_fill_msr.h"
  #include "include/dft_mesh.h"       
  #include "include/dft_fill_pde.h"
  #include "include/dft_mesh_lib.h"
  #include "include/dft_quadratue.h"
  #include "include/dft_fillp.h"
  #include "include/dft_mesh_surfaces.h"
  #include "include/dft_stencil.h"
  #include "include/dft_fill_rosen.h"
  #include "include/dft_newton.h"
  #include "include/dft_thermo.h"
  #include "include/dft_fill_rosen_rb.h"
  #include "include/dft_out_ads.h"
  #include "include/dft_vext.h"
  #include "include/dft_uww.h"
  #include "include/dft_images.h"
  #include "include/dft_potentials.h"
  #include "include/loca_con_bord.h"
  #include "include/loca_util.h"
#endif

