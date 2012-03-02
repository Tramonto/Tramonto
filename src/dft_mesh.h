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
extern void *ParameterList_list;
extern int Az_kspace;
extern int Max_gmres_iter;
#define SCREEN_BASIC       1
extern double Az_tolerance;
extern double Az_ilut_fill_param;
extern int Az_preconditioner;
extern int Az_scaling;
#define AM_taucs     22
#define AM_pardiso   21
#define AM_superludist 20
#define AM_superlu   19
#define AM_umfpack   18
#define AM_mumps     17
#define AM_klu       16
#define AM_lapack    15
extern int Az_solver;
void MY_read_update(int *N_update,int *update[],int N,int *nodes_x,int chunk,int input_option);
extern double *Lseg_IC;
extern int Nseg_IC;
extern double *Pore_rad_R_IC;
extern double X_const_mu;
#define OPTION_VARY  2
extern double *Pore_rad_L_IC;
#define PI    3.141592653589793238462643383279502884197169399375
#define OPTION_CYL   1
extern int Geom_flag;
#define FLAG_PBELEC -777
#define PB_ZONE      3 
#define FLAG_BULK   -888
#define BULK_ZONE    2 
void print_Nodes_to_zone(int *node_to_zone,char *Nodes2Zone_Filename);
int ijk_box_to_node_box(int *ijk_box);
extern int Nzone;
int element_to_node(int ielement);
int el_box_to_el(int iel_box);
void node_box_to_ijk_box(int node_box,int *ijk_box);
int element_box_to_node_box(int iel_box);
#define SMEAR_CHARGE 1
#define BACKGROUND   2
extern int Charge_type_local;
void els_charge_spheres(double radius,double *x,int *nelems,int *nelems_unique,int *elems,int charge_type);
#define NWALL_MAX_TYPE 50 
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
#define POINT_CHARGE 0
extern int Charge_type_atoms;
extern double *Charge_Diam;
extern double *Charge;
extern double **Charge_x;
#define NWALL_MAX 600 
extern double Elec_param_w[NWALL_MAX];
void print_charge_surf(double **charge_w_sum,char *output_file);
void print_charge_vol(double *charge_els,char *output_file);
void setup_volume_charge2(void);
void setup_linear_grad_of_charge(void);
void setup_volume_charge1(int iwall);
void bc_setup_const_charge(int iwall,int loc_inode);
extern int Nlocal_charge;
#define ATOMIC_CHARGE    3
#define CONST_CHARGE     2
extern int Type_bc_elec[NWALL_MAX_TYPE];
#define LAST_NODE_RESTART    4
#define LAST_NODE            3
double gsum_double(double c);
#define DOWN_FRONT   7
#define DOWN_BACK    3
#define LEFT_FRONT   6
#define LEFT_BACK    2
#define RIGHT_FRONT  5
#define RIGHT_BACK   1
#define UP_FRONT     4
#define UP_BACK      0
void surf_el_to_list(int loc_inode,int ilist,int *iel_box,int el,int type,int normal,int idim,double esize1,double esize2);
void find_local_els(int inode,int *iel,int *iel_box,int flag);
extern int Nnodes_per_el_S;
extern int Link[NWALL_MAX];
#define NDIM_MAX  3
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int **Link_list;
extern int *Nwall_this_link;
extern int **Xtest_reflect_TF;
int node_to_elem(int inode_all,int local_node,int *reflect_flag);
extern int ****Touch_domain_boundary;
#define WALL_EL        0 
#define FLUID_EL       1 
int el_to_el_box(int iel);
int node_to_elem_return_dim(int inode_all,int local_node,int *reflect_flag,int *idim_return,int *iside,int *periodic_flag);
extern int Nnodes_per_el_V;
#define atomic_centers                  3
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
void setup_zeroTF_and_Node2bound(FILE *fpecho,int ***el_type);
extern double **Charge_w_sum_els;
extern double *Charge_vol_els;
extern double **S_area_tot;
extern double ***S_area;
extern int ***Surf_elem_to_wall;
extern int **Surf_elem_type;
extern int ***Surf_normal;
extern int **Nelems_S;
void boundary_free(void);
void setup_surface_charge(FILE *fpecho);
extern int Surf_charge_flag;
extern int Vol_charge_flag;
extern int Ipot_wf_c;
void boundary_properties(FILE *fpecho);
void boundary_setup(char *file_echoinput);
int node_to_node_box(int inode);
#define MATRIX_FILL_NODAL 1   /* set to zero for physics based ordering */
int unk_box_to_unk(int i_box);
int node_box_to_node(int inode_box);
extern int Non_unique_G2B[4];
extern int Nunknowns_box;
extern int Elements_plane_box;
extern int Nodes_plane_box;
extern int Elements_x_box[NDIM_MAX];
extern int Nodes_x_box[NDIM_MAX];
extern int Pflag[3];
extern int Max_IJK[3];
extern int Max_IJK_box[3];
extern int Max_sten_length[3];
extern int Min_IJK[3];
extern int Min_IJK_box[3];
extern double X_1D_bc;
#define NCOMP_MAX 5
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
#define HARD_SPHERE  1
extern int Nunknowns;
extern double Area_surf_el[3];
extern double Vol_el;
extern int Elements_plane;
extern int Nodes_plane;
extern int Elements_x[NDIM_MAX];
extern int Nodes_x[NDIM_MAX];
extern double Esize_x[NDIM_MAX];
extern double Size_x[NDIM_MAX];
int round_to_int(double x);
extern int Nelements;
void setup_area_IC(void);
extern double *Area_IC;
#define DIFFUSIVE_INTERFACE 1
extern int Type_interface;
extern int Nnodes;
int gsum_int(int c);
extern int Coarser_jac;
void setup_wall_wall_potentials();
extern int Lprint_pmf;
void setup_external_field_n(int **nelems_w_per_w,int ***elems_w_per_w);
void read_external_field_n();
#define READ_VEXT_FALSE      0
extern int Restart_Vext;
void setup_zeroTF_and_Node2bound_new(FILE *fpecho,int ***el_type);
extern int *List_coarse_nodes;
extern int Nnodes_coarse_loc;
void set_mesh_coarsen_flag(void);
void zones_el_to_nodes(int *elem_zones);
extern int Nwall_Images;
extern int Imain_loop;
void setup_surface(FILE *fp2,int *nelems_f,int **nelems_w_per_w,int **elems_f,int ***elems_w_per_w,int *elem_zones,int ***el_type);
void setup_surface(FILE *fp2,int *nelems_f,int **nelems_w_per_w,int **elems_f,int ***elems_w_per_w,int *elem_zones,int ***el_type);
extern int Nnodes_wall_box;
extern int **List_wall_node;
extern int **Wall_touch_node;
extern int *Nwall_touch_node;
extern int *Nodes_wall_box;
extern int *Index_wall_nodes;
extern int Nlink;
extern double Dielec_bulk;
#define PERIODIC             1
#define REFLECT              2
extern int Type_bc[NDIM_MAX][2];
extern int Nelements_box;
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern int Ncomp;
#define IDEAL_GAS    0
extern int Ipot_ff_n;
extern int Lhard_surf;
extern int Nlists_HW;
extern int Nwall;
void linsolver_setup_control();
void node_to_position(int inode,double *NodePos);
void setup_basic_box(FILE *fpecho,int *update);
extern int *L2G_node;
extern int *B2L_node;
extern int *L2B_node;
extern int *B2G_unk;
extern double **Uww_link;
extern double **Uww;
extern double ***Xwall_delDOWN;
extern double ***Xwall_delUP;
extern double **X_wall;
#define VEXT_DIST_TO_SURF          2  /* any potential that is a function of only distance from the surface */
extern int Ipot_wf_n[NWALL_MAX_TYPE];
extern int Nwall_type;
extern int **Zero_density_TF;
extern double ***Vext_dash;
extern int Lvext_dash;
extern double **Vext;
extern double *Dielec;
#define COULOMB      1
extern int Ipot_ff_c;
extern int ***Wall_owners;
extern int **Nwall_owners;
extern int **Wall_elems;
extern int **Nodes_2_boundary_wall;
#define VERBOSE      3 
void safe_free(void **ptr);
void safe_free(void **ptr);
#if defined(DEC_ALPHA)
#define POW_INT powii
#endif
#if !(defined(DEC_ALPHA))
#define POW_INT (int)pow
#endif
extern int *Nodes_to_zone;
extern int Grad_dim;
#define FLAG_1DBC   -999
extern int *Mesh_coarsen_flag;
extern int Ndim;
void node_to_ijk(int node,int *ijk);
extern int *B2G_node;
extern int Nnodes_box;
extern void *LinProbMgr_manager;
void free_mesh_arrays(void);
void control_mesh(FILE *fpecho,char *output_file2,int print_flag,int *update);
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Mesh_coarsening;
#define LB_TIMINGS   3
#define LB_WEIGHTS   2
extern int Load_Bal_Flag;
extern int L1D_bc;
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int *Comm_offset_unk;
extern int *Comm_offset_node;
extern int *Comm_unk_proc;
extern int Num_Proc;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern int *Comm_node_proc;
void load_balance(int flag,double *fill_time,int *N_update,int **update);
extern int Nunk_per_node;
extern int Nnodes_per_proc;
void initialize_Aztec(int *N_update,int *update[]);
void setup_basic_domain(FILE *fpecho);
#define FILES_DEBUG        2
extern int Iwrite_files;
#define SCREEN_VERBOSE     3 
#define SCREEN_ERRORS_ONLY  0 
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void set_up_mesh(char *file_echoinput,char *output_file2);
