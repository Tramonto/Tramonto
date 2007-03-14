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
 *  File:  dft_globals.h
 *
 *  This file contains the declarations of some of the key global 
 *  varibles that are used often thoughout the code. This file should
 *  be included only once, in the main program file, and the related
 *  file dft_globals_const.h should be included in all other files that
 *  use one of these global variables.
 */

#include "dft_globals_const.h"

/* Basic Equation info */
int Phys2Nunk[NEQ_TYPE];  /* Number of unknowns of a particular equation type */
int Phys2Unk_first[NEQ_TYPE]; /* starting unknown number for a given equation type */
int Phys2Unk_last[NEQ_TYPE]; /* ending unknown number for a given equation type */
int Unk2Phys[3*NCOMP_MAX+NMER_MAX+NMER_MAX*NMER_MAX+13]; /* array that gives equation type
                                                                         given an unknown index */

/* Mesh info */

/*************** Global Mesh ********************************/
int     Ndim;            /* # of spatial dimensions of the problem      */
int     Nnodes;          /* # of nodes in the mesh                */
int     Nunknowns;       /* # of unknowns in the problem          */
int     Nelements;       /* # of elements in the mesh             */
int     Elements_plane;  /* # of elements in the x1-x2 plane      */
int     Nodes_plane;     /* # of nodes in the x1-x2 plane      */
int Nodes_x[NDIM_MAX];   /* Array[Ndim]: # nodes in idim       */
int Elements_x[NDIM_MAX];/* Array[Ndim]: # elements in idim    */
int     Max_sten_length[3];  /* The number of nodes in the longest stencil */

/************** Reference Variables *************************/
double Length_ref;
double Density_ref;
double Dielec_ref;
double Potential_ref;
double VEXT_MAX;

/*************** Extended Local Mesh ************************/
int     Nnodes_box;          /* # of nodes in the extended local mesh   */
int     Nunknowns_box;       /* # of unknowns in the extended local mesh  */
int     Nelements_box;       /* # of elements in the extended local mesh */
int     Elements_plane_box;  /* # of el in x1-x2 plane  of extended local mesh  */
int     Nodes_plane_box;     /* # of nodes in x1-x2 plane of extended local mesh*/
int Nodes_x_box[NDIM_MAX];   /* Array[Ndim]: # nodes in idim on extended local mesh */
int Elements_x_box[NDIM_MAX];/* Array[Ndim]: # elements in idim on extended local mesh */
int     Min_IJK_box[3]; /* The minimum IJK of the extended internals on this procesor */
int     Max_IJK_box[3];  /* The maximum IJK of the extended internals on this procesor */
int     Pflag[3];   /* PERIODIC flag array[Ndim]: TRUE Nodes_x_box=Nodes_x */

/************** Local Mesh **********************************/
int     Min_IJK[3];      /* The minimum IJK of the internals on this procesor */
int     Max_IJK[3];      /* The maximum IJK of the internals on this procesor */

/************** Communications arrays for Gather vectors ******/
int *Comm_node_proc;   /* array on proc 0 of nodes per processors */
int *Comm_unk_proc;   /* array on proc 0 of unknowns per processors */
int *Comm_offset_node; /* array on proc 0 of offsets of nodes  MPI_Gatherv*/
int *Comm_offset_unk;  /* array on proc 0 of offsets of unknowns MPI_Gatherv*/

/************** Other ***************************************/
double  Vol_el;          /* Volume of one element of our regular mesh        */
double Vol_in_surfs[NCOMP_MAX]; /* volume in all the surfaces of each list */
double  Area_surf_el[3];  /* Area of surface element w/ normal in direction idim*/
int     Nlists_HW;       /* Number of lists needed if hard walls (for mixtures)*/
int     Nel_wall_count;  /* Number of elements in the 0.5sigma of wall:rhobars*/
int    **Nelems_per_wall;   /*  Array[Nwall] of # of elements in a given wall */
int     Nnodes_per_proc;  /* Number of nodes owned by this processor         */
int     Nunk_int_and_ext; /* Number of unknownsneeded on this processor      */
int     *B2L_1stencil;     /* Box to Local array for nodes where rho_bar is needed*/
int     *B2L_unknowns;     /* Box to Local array for all unknowns */
int     *B2G_node;         /* Box to global array for all box nodes */
int     *B2G_unk;          /* Box to global array for all box unknowns */
int     *L2B_node;         /* Local to box  array for all local nodes */
int     *B2L_node;         /* Box to local  array for all local nodes */
int     *L2G_node;         /* Local to global coordinates  */
int     Nnodes_1stencil;  /* Number of nodes on proc where rho_bar is needed */
int     Nunk_per_node;   /* Number of unknowns per node (usually Ncomp       */
int     Nrho_bar;        /* Number of rhobar equations per node */
int     Nrho_bar_s;      /* Number of scalar rhobar equations per node */
int     Npoisson;      /* Number of scalar rhobar equations per node */
int     Ndiffusion;      /* Number of scalar rhobar equations per node */
int     NEL_unk;         /* Number of unknowns for the Euler-Lagrange equation */
int     Ntype_unk;       /* Number of equations defining segment types for polymer TC cases */
int     Nrho_bar_cavity; /* Number of nonlocal densities for the cavity function - WTC polymers */
int     Nrho_bar_bond; /* Number of nonlocal densities for the bond functionals - WTC polymers */
int     Nnodes_per_el_V;   /* Number of nodes per volume element             */
int     Nnodes_per_el_S;   /* Number of nodes per surface element            */
int     Plane_new_nodes; /* Indicates in which plane (xy,yz,xz)nodes are added*/
int     Pos_new_nodes;   /* Indicates where nodes are added                  */
double  Size_x[NDIM_MAX];    /* Array of the size of the domain in each dim. */
double  Esize_x[NDIM_MAX];   /* Array of the size of an element in each dim. */
int     Lmesh_refine;        /* Switch for automated mesh refinement         */
int     Type_bc[NDIM_MAX][2];/* Array of boundary conditions in each dim.    */
int     Non_unique_G2B[4];   /* Flag indicating which box dimensions 
                                wrap around completely, resulting in a
                                non-unique G2B mapping */
int   **Nodes_2_boundary_wall;  /* Array[Nlists_HW][Nnodes] -1 if nod b.n. else b.n.*/
int   **Wall_elems;  /* Array[Nlists_HW][Nelements] TRUE for wall elements  */
int   ****Touch_domain_boundary; /*Array[Nwall][Nlists_HW][[Ndim][2] =0 if surface hits 
                        left boundary = 1 if hits right domain boundary else -1 */
int   **Xtest_reflect_TF; /* Array[Nwall][Ndim] for reflections/wall-wall boundaries */ 
/* these are used to set up the Nodes_2_boundary_wall array then discarded*/
int   Nnodes_wall_box; /* Count number of nodes in box that touch a wall */
int  *Nodes_wall_box; /* Array to store which nodes touch a wall */
int  *Nwall_touch_node; /* Array to store number of walls touching a given node */
int  **Wall_touch_node; /*Array to store which walls touch a given node */
int  **List_wall_node; /*Array to store which walls touch a given node */
int  *Index_wall_nodes;  /*Array to store indexing in these mesh arrays*/

int First_time; /* for MSR preprocessing */

int     Nzone;          /* Number of diff. quadrature zones on the mesh      */
int    *Nodes_to_zone;   /* Array[Nnodes] of quadrature zones */
double  Rmax_zone[5];    /* Array distances from surfaces in quadrature zones */
int     Mesh_coarsening;  /* Flag indicating whether mesh coarsening is on */
int     Nnodes_coarse_loc; /* Number of coarse nodes local to a processor */
int    *List_coarse_nodes; /* List of coarse nodes local to a processor */
int    *Mesh_coarsen_flag;/* Flag (Nnodes) telling how much coarsening for a
                             given node, or negative values telling which
                             dimension to average over */
int     Coarser_jac;      /* Flag to switch on coarser jacobian than residual*/
double  Jac_grid;      /* Flag to switch on coarser jacobian than residual*/
int     Lcut_jac;  /* Logical to indicate if Jacobian stencils will be cut off */
double  Jac_threshold; /* Threshold level for Jacobian stencils ... max/Jac_threshold */

/* Continuation info */
int     Nodes_old;         /* # of nodes in the mesh of the previous run */
int Nodes_x_old[NDIM_MAX]; /* Array[Ndim]: # nodes in idim of previous run */
double *X_old;             /* Array of unknowns from previous run */
double *X2_old;             /* Array of unknowns from previous run */

int     Print_rho_type;  /* flag for printing file type of densities */
int     Print_rho_switch;  /* flag for printing of densities -- format*/
int     Lprint_gofr; /* flag for printing radial distribution functions */
int     Lprint_pmf; /* flag for printing radial distribution functions */
int     Print_mesh_switch;  /* flag for printing of densities -- format*/
int     Lper_area;  /*logical for per unit are outputs of params */ 
int     Lcount_reflect;  /*logical for per unit are outputs of params */ 
int     Nruns;           /* Number of runs to perform (varying the mesh)     */
double  Del_1[NWALL_MAX_TYPE];    /*Stepping parameter for field #1  */
double	Rho_max;         /* max rho when using an old solution for mesh contin*/
int     Imain_loop;   /* counter on the main loop of the program */
double  Scale_fac; /* for continuation in arrays !*/

/* Surface Physics info */
int     Nwall;           /* Number of surfaces in the calculation            */
int     Nwall_type;      /* Number of types of surface in the problem */
int     WallType[NWALL_MAX]; /*array containing type of each surface */
int     Nlink;           /* Number of macro-surfaces */
int     Link[NWALL_MAX]; /* Index iwall to a particular macro-surface */
int     *Nwall_this_link; /* Number of walls linked together in a particular macro-surface */
int     **Link_list; /* List of walls linked in a particular macro-surface */
int     Surface_type[NWALL_MAX_TYPE];    /* Type of surfaces of interest                     */
int     Orientation[NWALL_MAX_TYPE];  /* Orientation of planar/bumpy infinite walls*/
double  WallParam[NWALL_MAX_TYPE];/* Array[Nwall] of a characteristic wall parameter*/
double  WallParam_2[NWALL_MAX_TYPE];/* Array[Nwall] of a characteristic wall parameter*/
double  WallParam_3[NWALL_MAX_TYPE];/* Array[Nwall] of a characteristic wall parameter*/
double  WallParam_4[NWALL_MAX_TYPE];/* Array[Nwall] of a characteristic wall parameter*/
int     Lrough_surf[NWALL_MAX_TYPE]; /*Logical for rough surfaces */
double  Rough_precalc[NWALL_MAX_TYPE][MAX_ROUGH_BLOCK][MAX_ROUGH_BLOCK];
double  Rough_length[NWALL_MAX_TYPE];
double  WallPos[NDIM_MAX][NWALL_MAX]; /* Array of the centers of the surfaces*/

/* Fluid Physics info */
int     Ncomp;           /* Number of components in the current problem      */
double  Temp;            /* Reduced Temperature for LJ interactions          */
double  Temp_elec;       /* Reduced Temperature for Poisson's equation       */
int     Type_dielec;      /* choose how to handle dielectric constants in system */
double  Dielec_bulk;     /* the dielectric constant in the bulk fluid */
double  Dielec_pore;     /* the dielectric constant in the "pore" fluid */
double  Dielec_X;        /* distance from a surface defines "pore" fluid */
double *Dielec_wall;     /* the dielectric constant as function of wall type */
double *Dielec;          /* the dielectric constant as function of element box*/
int     Lpolarize[NCOMP_MAX]; /* Logical for if a species is polarizeable */
double  Pol[NCOMP_MAX];    /* Bulk polarization of each species */
double  Energy; /* surface free energy to return to Towhee  */


double  charge_fluid;    /* The charge in the fluid ... post processing */
double  Betap;           /* Pressure in units of kT sigma_ff[1]^3           */
double Betap_LBB;        /* Pressure for LBB of domain */
double Betap_RTF;        /* Pressure for RTF of domain */
double  Betap_id;        /* Ideal Gas Presseure in units of kT sigma_ff[1]^3 */
double  Betap_att;       /* Attractive Presseure in units of kT sigma_ff[1]^3 */
double  P_over_po;
double  Hs_diam[NCOMP_MAX]; /* Array of effective hard sphere diameters      */
int     L_isotherm;       /* Logical for isother vs. force vs. h data */
double  Rho_b[NCOMP_MAX];   /* Array[Ncomp] of component bulk values         */
double  Rhobar_b[10];   /* Array[Nrho_bar] of bulk rhobars        */
double  Rhobar_b_LBB[10];   /* Array[Nrho_bar] of bulk rhobars LBB    */
double  Rhobar_b_RTF[10];   /* Array[Nrho_bar] of bulk rhobars RTF    */
double  Dphi_Drhobar_b[10]; /* Array[Nrho_bar] of bulk energy derivs w/r/to rhobars      */
double  Dphi_Drhobar_LBB[10]; /* Array[Nrho_bar] of bulk energy derivs w/r/to rhobars LBB  */
double  Dphi_Drhobar_RTF[10]; /* Array[Nrho_bar] of bulk energy derivs w/r/to rhobars RTF  */
double  Rho_seg_b[NMER_MAX]; /* array of bulk segment densities */
double  Rho_seg_LBB[NMER_MAX];
double  Rho_seg_RTF[NMER_MAX];
double *Rhobar3_old;      /* Array[Nodes_box] of old values of rhobar3 */
double Xi_cav_b[4]; /* Array of bulk rhobars for cavity functions of WTC polymer functionals */
double Xi_cav_LBB[4]; /* Array of bulk rhobars for cavity functions of WTC polymer functionals */
double Xi_cav_RTF[4]; /* Array of bulk rhobars for cavity functions of WTC polymer functionals */
double BondWTC_b[NMER_MAX*NMER_MAX]; /* Array of bulk rhobars for bonds in WTC functionals*/
double BondWTC_LBB[NMER_MAX*NMER_MAX]; /* Array of bulk rhobars for bonds in WTC functionals*/
double BondWTC_RTF[NMER_MAX*NMER_MAX]; /* Array of bulk rhobars for bonds in WTC functionals*/

int     N_rho;   /* Number of density runs to perform         */
double  Del_rho;   /* How to step in density         */
double  Rho_coex[2];   /* Liquid-vapor coexisting densities         */
double  Betamu_hs_ex[NCOMP_MAX];/* Array of excess hardsphere chemical potentials*/
double  Betamu_ex_bondTC[NCOMP_MAX][NMER_MAX*NMER_MAX];/* Array of excess segment chemical potentials - WTC poolymer*/
double  Betamu_seg[NMER_MAX];/* Array of excess segment chemical potentials - WTC poolymer*/
double  Betamu_seg_LBB[NMER_MAX];/* Array of excess segment chemical potentials - WTC poolymer*/
double  Betamu_seg_RTF[NMER_MAX];/* Array of excess segment chemical potentials - WTC poolymer*/
double  Betamu[NCOMP_MAX]; /*Array[Ncomp] of chemical potentials */
double  Betamu_id[NCOMP_MAX]; /*Array[Ncomp] of ideal gas chemical potentials */
double  Betamu_wtc[NMER_MAX];
double  Betamu_wtc_LBB[NMER_MAX];
double  Betamu_wtc_RTF[NMER_MAX];

int     Ipot_ff_n;    /* Potential Type for neutral part of f-f interactions */
int     Ipot_wf_n[NWALL_MAX_TYPE];    /* Potential Type for neutral part of w-f interactions */
int     Type_pairPot;  /* Interaction potential to use for strict mean field DFT calculations*/
int     Type_vext1D;  /* Interaction potential to choose for external field calculations based on 1D potentials*/
int     Type_vext3D;  /* Interaction potential to choose for external field calculations based on 3D potentials*/
int     Type_uwwPot;  /* potential to use for computation of wall-wall interactions.  Used in 3D-atomic surface calculations */
int     Ipot_ww_n[NWALL_MAX_TYPE][NWALL_MAX_TYPE];    /* Potential Type for neutral part of w-f interactions */
int     Ipot_ff_c;    /* Potential Type for charged part of f-f interactions */
int     Ipot_wf_c;    /* Potential Type for charged part of w-f interactions */
int     Lhard_surf;   /* Logical indicating if surfaces are hard core */
int     Lvext_dash;   /* Logical indicating if the Vext_dash array should be set up */
int     Iliq_vap;     /* Type of Liquid-vapro profile                        */
int     Iguess1;        /* Type of initial guess */
int     Nsteps;         /* Number of steps for a step profile initial guess */
int     Orientation_step[NSTEPS_MAX]; /* orientation of the step profile */
double  Xstart_step[NSTEPS_MAX];  /* start position array for the step profile */
double  Xend_step[NSTEPS_MAX];  /* end position array for the step profile */
double  Rho_step[NCOMP_MAX][NSTEPS_MAX];  /* density array for a step profile */
int     Lbinodal;        /* Logical TF for binodal calculation */
double  Thickness;    /* Thickness parameter for wetting studies */
int     Mix_type;     /* Type of mixing rules */
double  Mass[NCOMP_MAX];           /* Array of the mass of each specie*/
double  Sigma_ff[NCOMP_MAX][NCOMP_MAX];/* Array of f-f interaction diameters */
double  Bond_ff[NCOMP_MAX][NCOMP_MAX];/* Array of f-f bond lengths for polymers */
double  Fac_overlap[NCOMP_MAX][NCOMP_MAX];/* Array of f-f bond lengths for polymers */
double  Fac_overlap_hs[NCOMP_MAX];/* Array of f-f bond lengths for polymers */
double  Eps_ff[NCOMP_MAX][NCOMP_MAX];  /* Array of f-f interaction energies  */
double  Cut_ff[NCOMP_MAX][NCOMP_MAX];  /* Array of f-f cutoff distances      */
double  Charge_f[NCOMP_MAX];           /* Array of the valence of each specie*/
double  Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];/* Array of w-f interaction diameters */
double  Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];  /* Array of w-f interaction energies  */
double  Cut_wf[NCOMP_MAX][NWALL_MAX_TYPE];  /* Array of w-f cutoff distances      */
double  Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];/* Array of w-w interaction diameters */
double  Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];  /* Array of w-w interaction energies  */
double  Cut_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];  /* Array of w-w cutoff distances      */
double  YukawaK_ff[NCOMP_MAX][NCOMP_MAX]; /* Yukawa decay constant on fluid-fluid interactions */
double  YukawaK_wf[NCOMP_MAX][NWALL_MAX_TYPE]; /* Yukawa decay constant on fluid-fluid interactions */
double  YukawaK_w[NWALL_MAX_TYPE]; /* Yukawa decay constant on fluid-fluid interactions */
double  YukawaK_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE]; /* Yukawa decay constant on fluid-fluid interactions */


int     **Lsemiperm;  /* Array of logicals for semipermeable surfaces */
double  **Vext_membrane; /*Array potentials for semi-perm surfaces */
double  **Vext_set;      /*Array of maximum set points for ext potentials */
double  **Vext;       /* External field array [Nnodes][Ncomp]           */
double  **Vext_coul;       /* External Coulomb field array [Nnodes]           */
double  *Vext_old;       /* For post processing External field array        */
double  ***Vext_dash; /* Derivative of external field [Nnodes][Ncomp][Nwall]  */
double  **Uww; /* wall-wall interactions [Nwall-1][Nwall-1]  */
double  **Uww_link; /* wall-wall interactions [Nlink-1][Nlink-1]  */
double  **X_wall; /* Distance from inode to iwall [Nnodes][Nwall]    */
double  **X_wall2; /* Distance from inode to iwall [Nnodes][Nwall]    */
int    **Zero_density_TF; /* array [Nnodes][icomp] indicates where VEXT_MAX */
double  Betamu_att[NCOMP_MAX];  /* sum over jcomp of Van der waals constant a(icomp,jcomp)*/
double  Avdw[NCOMP_MAX][NCOMP_MAX];    /* Van der waals constant a(icomp,jcomp)*/

double  Sigma_w[NWALL_MAX_TYPE];  /* Array[Nwall] of w-w interaction diameters    */
double  Eps_w[NWALL_MAX_TYPE];    /* Array[Nwall] of w-w interaction energies     */
double  Rho_w[NWALL_MAX_TYPE];    /* Array[Nwall] of w-w interaction energies     */
double  Elec_param_w[NWALL_MAX]; /* Array: surf charge (potential) per area  */
int     Type_bc_elec[NWALL_MAX_TYPE];/* Array of surface b.c.'s for charged systems*/
int     Nlocal_charge; /*Number of localized surface charges in the system */
double  **Charge_x;    /*Position of a given local charge [Nlocal_charge][Ndim]*/
double  *Charge_Diam; /*Diameter of a given local charge [Nlocal_charge]*/
double  *Charge; /*Value of the local charge [Nlocal_charge]*/
int     Charge_type_atoms; /* Type of charge distribution on charged atoms */
int     Charge_type_local; /* Type of charge distribution on added local charges */
double  *Deltac_b;   /* Array [icomp] of electrostatic correlations in bulk*/
double  **Charge_w_sum_els; /*Array[Nnodes_b][Ndim] of surface charge per area*/
double  *Charge_vol_els; /*Array[Nelements_box] of volumetric charge per elem*/
int     Vol_charge_flag; /* Flag for volumetric charges */
int     Surf_charge_flag; /* Flag for volumetric charges */
int     **Nelems_S; /* Array[Nlists_HW][Nodes_per_proc] of # surf elems the b.node touches*/
int     ***Surf_normal; /*Array[Nlists][Nodes_per_proc][Nelems_S] of unit normal vectors */
int     ***Surf_elem_to_wall; /*Array of wall to which a given surface element belongs...
                                 a given node may belong to more than one wall !! */
int     **Surf_elem_type; /*Array[Nodes_per_proc][Nelems_S] of surf elem type */
double  ***S_area; /*Array[Nlists][Nwall][Ndim] of total surf. area in idim on iwall*/
double  **S_area_tot; /*Array[Nlists][Nwall] total surf. area on iwall*/
int     **Nwall_owners; /*Array[Nilists][el_box] of number of walls 
                         (including images) that own a given element */
int     ***Wall_owners; /*Array[ilist][iel_box][Nwall_owners] that stores
                          all of the wall owners of a given element */

/* Steady State Solutions Info */

int    Lsteady_state;          /*True-False Steady State or Equilibrium Run*/
int    Linear_transport;       /*True-False Steady State or Equilibrium Run*/
double Velocity;               /* Convective velocity over diffusion coefficient*/
int    Grad_dim;               /*direction where gradient is implemented*/
int    L1D_bc;               /*logical for 1D boundary condition in Grad_dim */
double X_1D_bc;             /* distance where 1D boundary is applied */
double X_const_mu;             
double Rho_b_LBB[NCOMP_MAX];   /*Rho_b boundary condition left-bottom-back*/
double Rho_b_RTF[NCOMP_MAX];   /*Rho_b boundary condition right-top-front*/
double Elec_pot_LBB;           /*Electric potential boundary condition LBB*/
double Elec_pot_RTF;           /*Electric potential boundary condition RTF */
double Betamu_LBB[NCOMP_MAX];  /*Chemical Potential Boundary Condition LBB */
double Betamu_RTF[NCOMP_MAX];  /*Chemical Potential Boudary Condition RTF*/
double D_coef[NCOMP_MAX];  /*Diffusion Coefficients for ion species */
double *Pore_rad_L_IC;   /* array of left  Radii of ion chan pore segments (1D) */
double *Pore_rad_R_IC;   /* array of right Radii of ion chan pore segments(1D) */
double *Lseg_IC;         /* array of length of ion chan pore segments (1D)  */
int     Nseg_IC;         /* number of pore segments in an ion chan model */
double    *Area_IC;      /* 1D ion channel area per node (box units)*/
int   Geom_flag;   /* geometry flag for ion chann. see OPTION_ definitions*/

/* OUTPUT INTEGRAL PARAMETERS */
int    **Nel_hit;      /* number of elements hit by a given node in a given list */
int    **Nel_hit2;     /* same as prev. for a bulk fluid */
int    List[2];       /* which list numbers we care about for integrals*/
int    Imax;          /* how many lists are relevent to the case at hand 1 or 2 */
double Area;
double Fac_vol;
double Fac_area;
 

/* SOME CONSTANTS */
double Inv_4pi;               /* Precalculated value of 1/(4*pi)                    */
double Inv_rad[NCOMP_MAX];    /* Precalculated inverse of component radius          */
double Inv_4pir[NCOMP_MAX];   /* Precalculated inverse of component's 4 pi Radius   */
double Inv_4pirsq[NCOMP_MAX]; /* Precalculated inverse of component's 4 pi Radius^2 */

int     Type_func;    /* Type of functional for the calculation              */
int     Type_attr;    /* Type for handling attractions                       */
int     Type_coul;    /* Type for handling Coulombics                        */
int     Type_poly;    /* Type for handling polymers                          */
int     Lseg_densities; /* Logical to indicate that segement (rather than component) densities are treated in the code */
int     L_HSperturbation; /* Logical to indicate whether the run is base on perturbation of hard spheres */



/* Hard core type */
double   HS_diam[NCOMP_MAX];  /* Hard sphere diameters for the calculation */

/* Startup Info */
int     Restart;      /* Logical that switches between new prof & restart file*/
int     Restart_field[NEQ_TYPE];
int     Iprofile;     /* Specifies Liq-Solid,Vap-Solid,or Liq-Vap profile    */
double  Toler;        /* Tolerance for Newton-Rhapson iterations             */
int     Iwrite;       /* Do we want a complete or modified set of output data*/

/* Trilinos info */
void * LinProbMgr_manager;

/* Parallel Info, Aztec info */
int     Num_Proc; /* The total number of processors used in this calculation */
int     Proc;     /* The unique  processor number (from 0 to Num_Proc-1)     */
struct  Aztec_Struct Aztec; /* Structure to hold all the Aztec info          */
int     Load_Bal_Flag; /* Flag specifying type of laod balancing to do       */
int L_Schur; /* Switch to turn on Schur solvers */


/* Nonlinear Solver info */
int Max_Newton_iter;    /* Maximum # of Newton iterations (10 - 30)          */
double Newton_abs_tol,Newton_rel_tol; /* Convergence tolerances (update_soln)*/
double Min_update_frac; /* Minimum fraction to update solution to slow down
                           Newton's method */

/* Timers */
double Time_linsolver_first;
double Time_linsolver_av;
double Time_manager_first;
double Time_manager_av;
double Time_fill_first;
double Time_fill_av;


/* Linear Solver info */
int Az_solver;
int Az_kspace;
int Az_scaling;
int Az_preconditioner;
double Az_ilut_fill_param;
double Az_tolerance;
int Max_gmres_iter;

struct  Loca_Struct Loca; /* Information for continuation library */

/*
 * The global variable Stencil is a 3D array of Stencil_Struct
 * (for hard spheres of size [Nsten][Ncomp][Ncomp]). This variable is
 * defined in dft_stencil.h, and the extern statement is found here.
 */

struct Stencil_Struct ***Stencil;
int MPsten_Npts_R[NZONE_MAX];  /* # radial gauss pts. in MIDPOINT rule */
int MPsten_Npts_arc[NZONE_MAX]; /* # theta gauss pts. in MIDPOINT rule */
int MPsten_Npts_phi[NZONE_MAX]; /* # phi gauss pts. in MIDPOINT rule */
int Sten_Type[NSTEN]; /* on/off Flag for stencil types  */
int Sten_Choice_S[NSTEN][NZONE_MAX]; /* Quadrature type fore each stencil   */
int Sten_Choice_R[NSTEN][NZONE_MAX]; /* # Radial Gauss points for THETA_FNCs */
/****************************************************************************/
double Ads[NCOMP_MAX][2];
double Ads_ex[NCOMP_MAX][2];
double *Integration_profile; /* a place to put the integrand as a function of position */
/****************************************************************************/

/* Polymer variables */
double Deltar_cr,Gauss_a,Gauss_k,***Rism_cr;
double Crfac;
double Cr_rad[NCOMP_MAX][NCOMP_MAX];
double Cr_rad_hs[NCOMP_MAX][NCOMP_MAX];
double Bupdate_fact;
int Nblock[NCOMP_MAX],Ntype_mer,Nmer[NCOMP_MAX],Type_mer[NCOMP_MAX][NMER_MAX];
int Npol_comp,Nmer_t[NCOMP_MAX][NBLOCK_MAX],Last_nz_cr;
int Nmer_t_total[NBLOCK_MAX];
int Nseg_tot;
int Nseg_type[NCOMP_MAX];
int Bupdate_iters,Geqn_start[NCOMP_MAX];
char Cr_file[40],Cr_file2[40],Cr_file3[40],Cr_file4[40];
double Cr_break[2];
int  Ncr_files;
int *Unk_to_Poly, *Unk_to_Seg, *Unk_to_Bond, ***Poly_to_Unk, **Poly_to_Unk_SegAll;
int Ngeqn_tot, Nbonds, **Nbond,***Bonds; 
int *Pol_Sym;
int *Pol_Sym_Seg;
int *BondAll_to_isegAll;
int *BondAll_to_ibond;
int Unk2Comp[NMER_MAX],SegChain2SegAll[NCOMP_MAX][NMER_MAX],**Bonds_SegAll,*Nbonds_SegAll;

