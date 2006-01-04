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
 *  FILE: dft_stencil.c
 *
 *  This file contains the routines to precalculate the "stencils,"
 *  which are the integration schemes for the non-local interactions.
 *  The results will be strored in the global variable Stencil, which
 *  is a 2D array of Stencil_Struct (defined in dft_stencil_const.h).
 *
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"

/* Prototypes for functions found in this file */
static void merge_load_stencil(struct Stencil_Struct *, int **, double *,int *, int *);
static void print_out_stencil(int, int,int, int, FILE *);
static void sort_stencil(struct Stencil_Struct *);
static double calc_sten_rad(int , int , int);
static double calc_sten_vol(int , int , int);
static double get_weight_from_stencil(int , int , int , double , double,
                                      int , double *, double *);
static int calc_in_out_on(double *, double *, double );
static void renormalize_stencil(struct Stencil_Struct *, double);
double int_cr(double,double,double,int,int, int, double, double *);
double gauss(double, int, int);
int ijk_to_isten_index(int *,int *);
void shorten_stencil(struct Stencil_Struct *);



void calc_stencils(void)
{
 /*
  * Local variable declarations
  */

  char *yo = "calc_stencils";
  int izone, isten, 
      icomp, i, j, k;   /* counters in for loops         */
  int jcomp, imax, jmax, idim, ig, jg, kg, iii; /* counters in for loops      */

  int    **el_offsets;
  double *el_weights;

  int max_sten_nodes, el_in_radius[3], t1=0, in_out_on_flag;
  int nsten_max;
  double vol_sten, sten_rad, el_vol;

  int ix;
  double sum_w;
 
  /* quadrature and interpolation stuff */
  double weight, x_left[3], x_right[3];
  double qp[3], point[3], ok_distance=0.0;
  int ngp, ngp1, ngp2, ngp3, npt=0;
  double gp1[12],gp2[12],gp3[12],gw1[12],gw2[12],gw3[12];
  double *gp, *gw;
  double r_center_sq,x_mid[3];
  double radius,radius_sq, inv_sten_rad_sq, inv_npt;


  /* u_attract quadrature off-dimensions*/
  int ngpu, ngpu1,ngpu2,ngpu3;
  double gpu1[40], gwu1[40], gpu2[40], gwu2[40], gpu3[40],gwu3[40];
  double *gpu, *gwu;

  struct Stencil_Struct *sten;
  int *index_sten;
  FILE *ifp=NULL, *ifp2=NULL;

  double esize_zone[3];
  int  zone_coarseness=0;

  /********************** BEGIN EXECUTION ************************************/

  if (Proc==0) {
    if (Iwrite!=NO_SCREEN){
         printf("\n---------------------------------------------------------------------------------\n");
         printf("%s: Calculating stencils ... \n",yo);
    }
    t1 = MPI_Wtime();
    if (Iwrite == VERBOSE) ifp = fopen("stencil.out", "w");
    ifp2 = fopen("dft_out.lis","a+");
  }

 /*
  * Allocate the Stencil variable to be a 3D array of structures
  * Now it can be accessed as Stencil[i][j][k].Length
  */

  Stencil = (struct Stencil_Struct ***) array_alloc
             (3, NSTEN, Nzone, Ncomp*Ncomp, sizeof(struct Stencil_Struct));

  el_offsets = (int **) array_alloc (2, Nnodes_per_el_V, Ndim, sizeof(int));
  el_weights = (double *) array_alloc (1, Nnodes_per_el_V, sizeof(double));

  ngp1  = 6; ngp2  = 3; ngp3 = 1;
  ngpu1 = 40; ngpu2 = 20; ngpu3 = 12;

  set_gauss_quad(ngp1, gp1, gw1);
  set_gauss_quad(ngp2, gp2, gw2);
  set_gauss_quad(ngp3, gp3, gw3);

  /* zero the maximum stencil length .... used later for
     setting up local coordinate systems */
  for (idim=0; idim<Ndim; idim++) {
       Max_sten_length[idim] = 0;
       Sten_length_hs[idim] = 0;
  }
 
 /* Loop over each quadrature zone */
  for (izone=0; izone<Nzone; izone++){

  if (Coarser_jac ==5 && izone == Nzone-1 ){
    for (idim=0; idim<Ndim; idim++){
        zone_coarseness = (int)(Jac_grid/Esize_x[idim]+0.000000001);
        esize_zone[idim]=Esize_x[idim]* zone_coarseness;
    }
  }
  else{
    zone_coarseness = POW_INT(2,izone);
    for (idim=0; idim<Ndim; idim++)
       esize_zone[idim]=Esize_x[idim]* zone_coarseness;
  }

 /* Loop over each stencil type, only proceed if stencil type is true */
 
  for (isten=0; isten<NSTEN; isten++)  
  if (Sten_Type[isten]) {

 
    if (isten == POLYMER_GAUSS) imax = 1;  /* currently a single Gauss bl must
                                              be used for all polymer segments*/
    else                   imax = Ncomp;

    /* here we split stencils into those that depend only on one species and those that
       depend on pairs of species */
    if ((isten == DELTA_FN && (Type_poly==NONE || Type_poly==WTC)) || isten==THETA_FN ||isten==THETA_FN_SIG) jmax=1;
    else if (isten == U_ATTRACT || isten == THETA_CHARGE || isten == POLYMER_CR ||
        (isten == DELTA_FN && Type_poly!=NONE && Type_poly != WTC) || isten==DELTA_FN_BOND) jmax = Ncomp;
    else {
       printf("problems defining the number of stencil functions to compute for the chosen type\n!");
       exit(-1);
    }

    if (isten == U_ATTRACT || isten == THETA_CHARGE || isten == POLYMER_CR){
       if (isten== THETA_CHARGE) ngpu = 20;
       set_gauss_quad(ngpu1, gpu1, gwu1);
       set_gauss_quad(ngpu2, gpu2, gwu2);
       set_gauss_quad(ngpu3, gpu3, gwu3);
    }

    for (icomp=0; icomp<imax; icomp++) {
    for (jcomp=0; jcomp<jmax; jcomp++) {

      /* calculate stencil volume (e.g. 4/3 pi r^3), radius, and element vol */

      vol_sten = calc_sten_vol(isten, icomp, jcomp);
      sten_rad = calc_sten_rad(isten, icomp, jcomp);
      inv_sten_rad_sq = 1.0/(sten_rad*sten_rad);
      el_vol=1.0; for (i=0; i<Ndim; i++) el_vol *= esize_zone[i];

      /*
       * Figure out elements that form a bounding box for the stencil 
       * and the maximum number of nodes in the stencil
       */

      max_sten_nodes = 1.0;
      for (i=0; i<Ndim; i++) {
        el_in_radius[i] = (int) ((float)sten_rad/(float) esize_zone[i] + 0.999);
        max_sten_nodes *= 2 * el_in_radius[i] + 1;
	
        if (isten == U_ATTRACT || isten == THETA_CHARGE 
            || isten == POLYMER_CR || isten == POLYMER_GAUSS) 
                            nsten_max = el_in_radius[i];
        else{               nsten_max = 2*el_in_radius[i];}
        if (isten == DELTA_FN || isten==THETA_FN){
                if (el_in_radius[i]*zone_coarseness > Sten_length_hs[i])
                    Sten_length_hs[i] = el_in_radius[i]*zone_coarseness;
        }
        if (nsten_max*zone_coarseness > Max_sten_length[i]) 
              Max_sten_length[i] = nsten_max*zone_coarseness;
      }

      /* Allocate current stencil using this upper bound length */

      sten = &(Stencil[isten][izone][icomp+Ncomp*jcomp]);
      sten->Length = 0;
      sten->Offset = (int **) array_alloc(2, max_sten_nodes, Ndim, sizeof(int));
      sten->Weight = (double *) array_alloc(1, max_sten_nodes, sizeof(double));
      if (Lhard_surf) sten->HW_Weight = (double **)
           array_alloc(2, max_sten_nodes, Nnodes_per_el_V, sizeof(double));

      /* Allocate temporary index array to keep track of which 
         nodes have been set and which have not. -1 is the flag value
         for not yet set ! */
      index_sten = (int *) array_alloc(1, max_sten_nodes, sizeof(int));
      for (i=0; i<max_sten_nodes; i++) index_sten[i] = -1;
        

      /* Loop through all elements in bounding box, and perform quadrature */

      switch (Ndim) {

        case 1:

          /* loop over all elements, which contain nodes i and i+1 */

          for (i = -el_in_radius[0]; i < el_in_radius[0]; i++) {
             x_left[0]  = i * esize_zone[0];
             x_right[0] = (i+1) * esize_zone[0];
             el_offsets[0][0] = i;
             el_offsets[1][0] = i + 1;
             x_mid[0] = x_left[0] + 0.5*esize_zone[0];
             r_center_sq = (x_mid[0]*x_mid[0]);
             if (r_center_sq <= 4.0){
                 if (isten == U_ATTRACT)
                     { ngp = ngp3; gp = &gp3[0]; gw = &gw3[0];}
                 else{ ngp = ngp1; gp = &gp1[0]; gw = &gw1[0];}
                 ngpu = ngpu1; gpu = &gpu1[0]; gwu = &gwu1[0];
             }
             else if (r_center_sq <= 16.0){
                 if (isten == U_ATTRACT)
                     { ngp = ngp3; gp = &gp3[0]; gw = &gw3[0];}
                 else{ ngp = ngp2; gp = &gp2[0]; gw = &gw2[0];}
                 ngpu = ngpu2; gpu = &gpu2[0]; gwu = &gwu2[0];
             }
             else{
                 if (isten == U_ATTRACT)
                     { ngp = ngp3; gp = &gp3[0]; gw = &gw3[0];}
                 else{ ngp = ngp2; gp = &gp2[0]; gw = &gw2[0];}
                 ngpu = ngpu3; gpu = &gpu3[0]; gwu = &gwu3[0];
             }

             for (iii=0; iii < Nnodes_per_el_V; iii++)
                el_weights[iii] = 0.0;

             /* test if element is out, in, or stradles radius */

             in_out_on_flag = calc_in_out_on(x_left, x_right, sten_rad);
               
             if (in_out_on_flag == 1) {
               /* nothing to do -- element outside of stencil radius */
             }
             else{
               if (in_out_on_flag == -1) {
                 /* element completely within stencil radius -- use Gauss Quad*/

                 for (ig=0; ig < ngp; ig++) {

                   point[0] = x_left[0] + gp[ig] * esize_zone[0];
  
                   /*radius = fabs(point[0]) / sten_rad;*/

                   radius_sq = point[0]*point[0]*inv_sten_rad_sq;

                   if (isten != POLYMER_CR)
                     weight = get_weight_from_stencil(isten, icomp, jcomp, 
                                          radius_sq,sten_rad, ngpu, gpu, gwu);
                   else {
                     radius = fabs(point[0]);
                     weight = get_weight_from_stencil(isten, icomp, jcomp, 
                                          radius,sten_rad, ngpu, gpu, gwu);
                   }
                   weight *= gw[ig] * el_vol;

                   el_weights[0] += weight * (1.0-gp[ig]); 
                   el_weights[1] += weight * gp[ig]; 
  
                 }
  
               } 
               else {
                 /* element stradles stencil boundary -- use lots of dumb quad*/

                 switch (isten) {
                    case DELTA_FN:      npt = 150; break; /* singular at r=1 */
                    case THETA_FN:      npt =  40; break;
                    case U_ATTRACT:     npt =  20; break;
                    case THETA_CHARGE:  npt =  20; break;
                    case POLYMER_CR:    npt =  20; break;
                    case POLYMER_GAUSS:      npt =  20; break;
                    case THETA_FN_SIG:  npt = 40; break;
                    case DELTA_FN_BOND: npt = 150; break;
                 }
                 inv_npt = 1.0 / (double) npt;
                 for (ig=0; ig < npt; ig++) {
                   
                   qp[0] =  (ig + 0.5) * inv_npt;

                   point[0] = x_left[0] + qp[0] * esize_zone[0];
  
                   /*radius = fabs(point[0]) / sten_rad;*/
                   radius_sq = point[0]*point[0]*inv_sten_rad_sq;

                   if (radius_sq <= 1.000001) {
                     if (isten != POLYMER_CR)
                       weight = get_weight_from_stencil(isten, icomp, jcomp, 
                                            radius_sq,sten_rad, ngpu, gpu, gwu);
                     else {
                       radius = fabs(point[0]);
                       weight = get_weight_from_stencil(isten, icomp, jcomp, 
                                            radius,sten_rad, ngpu, gpu, gwu);
                     }
                     if (fabs(radius_sq-1.0) <= 1.e-8) weight *= 0.5;
                     weight *= inv_npt * el_vol; 

                     el_weights[0] += weight * (1.0-qp[0]); 
                     el_weights[1] += weight * qp[0]; 
                   }
                 }
               
               }

               merge_load_stencil(sten, el_offsets, el_weights,el_in_radius,index_sten);
             }
          }
          break;

        case 2:

          /* loop over all elements, which contain nodes i and i+1 */

          for (i = -el_in_radius[0]; i < el_in_radius[0]; i++) {
             x_left[0]  = i * esize_zone[0];
             x_right[0] = (i+1) * esize_zone[0];
             el_offsets[0][0] = i;
             el_offsets[1][0] = i + 1;
             el_offsets[2][0] = i;
             el_offsets[3][0] = i + 1;
          for (j = -el_in_radius[1]; j < el_in_radius[1]; j++) {
             x_left[1]  = j * esize_zone[1];
             x_right[1] = (j+1) * esize_zone[1];
             el_offsets[0][1] = j;
             el_offsets[1][1] = j;
             el_offsets[2][1] = j + 1;
             el_offsets[3][1] = j + 1;

             x_mid[0] = x_left[0] + 0.5*esize_zone[0];
             x_mid[1] = x_left[1] + 0.5*esize_zone[1];

             r_center_sq = (x_mid[0]*x_mid[0]) + (x_mid[1]*x_mid[1]);
             if (r_center_sq <= 4.0000001){
                 if (isten == U_ATTRACT)
                     { ngp = ngp3; gp = &gp3[0]; gw = &gw3[0];}
                 else{ ngp = ngp1; gp = &gp1[0]; gw = &gw1[0];}
                 ngpu = ngpu1; gpu = &gpu1[0]; gwu = &gwu1[0];
             }
             else if (r_center_sq <= 16.0000001){
                 if (isten == U_ATTRACT)
                     { ngp = ngp3; gp = &gp3[0]; gw = &gw3[0];}
                 else{ ngp = ngp2; gp = &gp2[0]; gw = &gw2[0];}
                 ngpu = ngpu2; gpu = &gpu2[0]; gwu = &gwu2[0];
             }
             else{
                 if (isten == U_ATTRACT)
                     { ngp = ngp3; gp = &gp3[0]; gw = &gw3[0];}
                 else{ ngp = ngp2; gp = &gp2[0]; gw = &gw2[0];}
                 ngpu = ngpu3; gpu = &gpu3[0]; gwu = &gwu3[0];
             }

             for (iii=0; iii < Nnodes_per_el_V; iii++)
                el_weights[iii] = 0.0;

             /* test if element is out, in, or stradles radius */

             in_out_on_flag = calc_in_out_on(x_left, x_right, sten_rad);
               
             if (in_out_on_flag == 1) {
               /* nothing to do -- element outside of stencil radius */
             }
             else {
               if (in_out_on_flag == -1) {
                 /* element completely within stencil radius -- use Gauss Quad*/

                 for (ig=0; ig < ngp; ig++) {

                   point[0] = x_left[0] + gp[ig] * esize_zone[0];
  
                   for (jg=0; jg < ngp; jg++) {
                     point[1] = x_left[1] + gp[jg] * esize_zone[1];
                     /*radius = sqrt(point[0]*point[0] + point[1]*point[1])
                                                               / sten_rad;*/
                     radius_sq = (point[0]*point[0] + point[1]*point[1])
                                 *inv_sten_rad_sq;

                   if (isten != POLYMER_CR)
                     weight = get_weight_from_stencil(isten, icomp, jcomp, 
                                          radius_sq,sten_rad, ngpu, gpu, gwu);
                   else {
                     radius = sqrt(point[0]*point[0] + point[1]*point[1]);
                     weight = get_weight_from_stencil(isten, icomp, jcomp, 
                                          radius,sten_rad, ngpu, gpu, gwu);
                   }

                     weight *= gw[ig] * gw [jg] * el_vol;

                     el_weights[0] += weight * (1.0-gp[ig]) * (1.0-gp[jg]); 
                     el_weights[1] += weight * gp[ig] * (1.0-gp[jg]); 
                     el_weights[2] += weight * (1.0-gp[ig]) * gp[jg]; 
                     el_weights[3] += weight * gp[ig] * gp[jg]; 
  
                   }
                 }
  
               } 
               else {
                 /* element stradles stencil boundary -- use lots of dumb quad*/

                 switch (isten) {
                    case DELTA_FN:      npt = 150; break; /* singular at r=1 */
                    case THETA_FN:      npt =  40; break;
                    case U_ATTRACT:     npt =  20; break;
                    case THETA_CHARGE:  npt =  20; break;
                    case POLYMER_CR:    npt =  20; break;
                    case POLYMER_GAUSS:      npt =  20; break;
                 }
                 inv_npt = 1.0 / (double) npt;

                 for (ig=0; ig < npt; ig++) {
                   
                   qp[0] =  (ig + 0.5) * inv_npt;

                   point[0] = x_left[0] + qp[0] * esize_zone[0];
  
                   for (jg=0; jg < npt; jg++) {

                     qp[1] =  (jg + 0.5) * inv_npt;
                     point[1] = x_left[1] + qp[1] * esize_zone[1];
                     /*radius = sqrt(point[0]*point[0] + point[1]*point[1])
                                                               / sten_rad;*/
                     radius_sq = (point[0]*point[0] + point[1]*point[1])
                                 * inv_sten_rad_sq;

                     if (radius_sq <= 1.0000001) {

                   if (isten != POLYMER_CR)
                     weight = get_weight_from_stencil(isten, icomp, jcomp, 
                                          radius_sq,sten_rad, ngpu, gpu, gwu);
                   else {
                     radius = sqrt(point[0]*point[0] + point[1]*point[1]);
                     weight = get_weight_from_stencil(isten, icomp, jcomp, 
                                          radius,sten_rad, ngpu, gpu, gwu);
                   }
                       if (fabs(radius_sq-1.0) <= 1.e-8) weight *= 0.5;
                       weight *= inv_npt * inv_npt * el_vol; 
  
                       el_weights[0] += weight * (1.0-qp[0]) * (1.0-qp[1]); 
                       el_weights[1] += weight * qp[0] * (1.0-qp[1]); 
                       el_weights[2] += weight * (1.0-qp[0]) * qp[1]; 
                       el_weights[3] += weight * qp[0] * qp[1]; 
                     }
                   }
                 }
               
               }

               merge_load_stencil(sten, el_offsets, el_weights,el_in_radius,index_sten);
             }
          }
          }
          break;

        case 3:

          /* loop over all elements, which contain nodes i and i+1 */

          for (i = -el_in_radius[0]; i < el_in_radius[0]; i++) {
             x_left[0]  = i * esize_zone[0];
             x_right[0] = (i+1) * esize_zone[0];
             el_offsets[0][0] = i;
             el_offsets[1][0] = i + 1;
             el_offsets[2][0] = i;
             el_offsets[3][0] = i + 1;
             el_offsets[4][0] = i;
             el_offsets[5][0] = i + 1;
             el_offsets[6][0] = i;
             el_offsets[7][0] = i + 1;
          for (j = -el_in_radius[1]; j < el_in_radius[1]; j++) {
             x_left[1]  = j * esize_zone[1];
             x_right[1] = (j+1) * esize_zone[1];
             el_offsets[0][1] = j;
             el_offsets[1][1] = j;
             el_offsets[2][1] = j + 1;
             el_offsets[3][1] = j + 1;
             el_offsets[4][1] = j;
             el_offsets[5][1] = j;
             el_offsets[6][1] = j + 1;
             el_offsets[7][1] = j + 1;
          for (k = -el_in_radius[2]; k < el_in_radius[2]; k++) {
             x_left[2]  = k * esize_zone[2];
             x_right[2] = (k+1) * esize_zone[2];
             el_offsets[0][2] = k;
             el_offsets[1][2] = k;
             el_offsets[2][2] = k;
             el_offsets[3][2] = k;
             el_offsets[4][2] = k + 1;
             el_offsets[5][2] = k + 1;
             el_offsets[6][2] = k + 1;
             el_offsets[7][2] = k + 1;

             x_mid[0] = x_left[0] + 0.5*esize_zone[0];
             x_mid[1] = x_left[1] + 0.5*esize_zone[1];
             x_mid[2] = x_left[2] + 0.5*esize_zone[2];

             r_center_sq = (x_mid[0]*x_mid[0]) + (x_mid[1]*x_mid[1] + 
                            x_mid[2]*x_mid[2]);
             if (r_center_sq <= 4.0){
                 if (isten == U_ATTRACT)
                     { ngp = ngp3; gp = &gp3[0]; gw = &gw3[0];}
                 else{ ngp = ngp1; gp = &gp1[0]; gw = &gw1[0];}
                 ngpu = ngpu1; gpu = &gpu1[0]; gwu = &gwu1[0];
             }
             else if (r_center_sq <= 16.0){
                 if (isten == U_ATTRACT)
                     { ngp = ngp3; gp = &gp3[0]; gw = &gw3[0];}
                 else{ ngp = ngp2; gp = &gp2[0]; gw = &gw2[0];}
                 ngpu = ngpu2; gpu = &gpu2[0]; gwu = &gwu2[0];
             }
             else{
                 if (isten == U_ATTRACT)
                     { ngp = ngp3; gp = &gp3[0]; gw = &gw3[0];}
                 else{ ngp = ngp2; gp = &gp2[0]; gw = &gw2[0];}
                 ngpu = ngpu3; gpu = &gpu3[0]; gwu = &gwu3[0];
             }

             for (iii=0; iii < Nnodes_per_el_V; iii++)
                el_weights[iii] = 0.0;

             /* test if element is out, in, or stradles radius */

             in_out_on_flag = calc_in_out_on(x_left, x_right, sten_rad);
               
             if ( (in_out_on_flag == 1) ||
                  ((in_out_on_flag == -1) && (isten == DELTA_FN))){
               /* nothing to do in this case -- element outside of stencil   */
               /* radius or completely inside of a delta function stencil    */
             }
             else {
               if (in_out_on_flag == -1) {
                 /* element completely within stencil radius -- use Gauss Quad*/

                 for (ig=0; ig < ngp; ig++) {

                   point[0] = x_left[0] + gp[ig] * esize_zone[0];
  
                   for (jg=0; jg < ngp; jg++) {

                     point[1] = x_left[1] + gp[jg] * esize_zone[1];
  
                     for (kg=0; kg < ngp; kg++) {
  
                       point[2] = x_left[2] + gp[kg] * esize_zone[2];
                       /*radius = sqrt(point[0]*point[0] + point[1]*point[1] 
                                     + point[2]*point[2]) / sten_rad;*/
                       radius_sq = (point[0]*point[0] + point[1]*point[1]
                                     + point[2]*point[2])*inv_sten_rad_sq;
  
                   if (isten != POLYMER_CR)
                     weight = get_weight_from_stencil(isten, icomp, jcomp, 
                                          radius_sq,sten_rad, ngpu, gpu, gwu);
                   else {
                     radius = sqrt(point[0]*point[0] + point[1]*point[1]
                                                     + point[2]*point[2]);
                     weight = get_weight_from_stencil(isten, icomp, jcomp, 
                                          radius,sten_rad, ngpu, gpu, gwu);
                   }

                       weight *= gw[ig] * gw [jg] * gw[kg] * el_vol;
  
                       el_weights[0] += weight *
                                    (1.0-gp[ig]) * (1.0-gp[jg]) * (1.0-gp[kg]);
                       el_weights[1] += weight *
                                         gp[ig]  * (1.0-gp[jg]) * (1.0-gp[kg]);
                       el_weights[2] += weight *
                                    (1.0-gp[ig]) *      gp[jg]  * (1.0-gp[kg]);
                       el_weights[3] += weight *
                                         gp[ig] *       gp[jg]  * (1.0-gp[kg]); 
                       el_weights[4] += weight *
                                    (1.0-gp[ig]) * (1.0-gp[jg]) *      gp[kg];
                       el_weights[5] += weight *
                                         gp[ig]  * (1.0-gp[jg]) *      gp[kg];
                       el_weights[6] += weight *
                                    (1.0-gp[ig]) *      gp[jg]  *      gp[kg];
                       el_weights[7] += weight *
                                         gp[ig]  *      gp[jg]  *      gp[kg]; 
                     }
                   }
                 }
  
               } 
               else {
                 /* element stradles stencil boundary -- use lots of dumb quadrature */

                 switch (isten) {
                    case DELTA_FN:      npt =  40; break;
                    case THETA_FN:      npt =  20; break;
                    case U_ATTRACT:     
                        if (Ndim==3)    npt = ngp3; 
                        else            npt =  10; 
                        break;
                    case THETA_CHARGE:  npt =  20; break;
                    case POLYMER_CR:    npt =  20; break;
                    case POLYMER_GAUSS:      npt =  20; break;
                 }
                 inv_npt = 1.0 / (double) npt;

                 /* define distance from the r=1 that counts as being on sphere's surface */

                 if (isten == DELTA_FN) 
                    ok_distance = sqrt(esize_zone[0]*esize_zone[0] +
                                       esize_zone[1]*esize_zone[1] +
                                       esize_zone[2]*esize_zone[2])
                                        *(0.5 * inv_npt);

                 for (ig=0; ig < npt; ig++) {
                   
                   qp[0] =  (ig + 0.5) * inv_npt;

                   point[0] = x_left[0] + qp[0] * esize_zone[0];
  
                   for (jg=0; jg < npt; jg++) {

                     qp[1] =  (jg + 0.5) * inv_npt;
                     point[1] = x_left[1] + qp[1] * esize_zone[1];
  
                     for (kg=0; kg < npt; kg++) {

                       qp[2] =  (kg + 0.5) * inv_npt;
                       point[2] = x_left[2] + qp[2] * esize_zone[2];
                       /*radius = sqrt(point[0]*point[0] + point[1]*point[1]
                                     + point[2]*point[2]) / sten_rad;*/
                       radius_sq = (point[0]*point[0] + point[1]*point[1]
                                     + point[2]*point[2])*inv_sten_rad_sq;

                       if ( ((radius_sq <= 1.0) && (isten != DELTA_FN)) || 
                            ((isten == DELTA_FN) && (fabs(radius_sq-1.0) < 2*ok_distance)) ) {
                   if (isten != POLYMER_CR)
                     weight = get_weight_from_stencil(isten, icomp, jcomp, 
                                          radius_sq,sten_rad, ngpu, gpu, gwu);
                   else {
                     radius = sqrt(point[0]*point[0] + point[1]*point[1]
                                                     + point[2]*point[2]);
                     weight = get_weight_from_stencil(isten, icomp, jcomp, 
                                          radius,sten_rad, ngpu, gpu, gwu);
                   }

                         if (fabs(radius_sq-1.0) <= 1.e-8 && isten != DELTA_FN) weight *= 0.5;

                         weight *=  (inv_npt*inv_npt*inv_npt) * el_vol; 
    
                         el_weights[0] += weight * 
                                        (1.0-qp[0]) * (1.0-qp[1]) * (1.0-qp[2]);
                         el_weights[1] += weight *
                                             qp[0]  * (1.0-qp[1]) * (1.0-qp[2]);
                         el_weights[2] += weight *
                                        (1.0-qp[0]) *      qp[1]  * (1.0-qp[2]);
                         el_weights[3] += weight *
                                             qp[0]  *      qp[1]  * (1.0-qp[2]);
                         el_weights[4] += weight *
                                        (1.0-qp[0]) * (1.0-qp[1]) *      qp[2];
                         el_weights[5] += weight *
                                             qp[0]  * (1.0-qp[1]) *      qp[2];
                         el_weights[6] += weight *
                                        (1.0-qp[0]) *      qp[1] *       qp[2];
                         el_weights[7] += weight *
                                             qp[0]  *      qp[1] *       qp[2];
                       }
                     }
                   }
                 }
               
               }
               merge_load_stencil(sten, el_offsets, el_weights,el_in_radius,index_sten);
             }
          }
          }
          }
          break;

      } 

      safe_free((void *)&index_sten);

      /* Sort Stencil in ascending order by node number */
  
/*  SORT DOESN"T SORT HW STUFF -- DO NOT USE AS IS    
 * sort_stencil(sten);
 */

      /* coarse stencils for big esize_zone must be scaled to real Esize*/

      if (izone > 0){
        for (i=0; i<sten->Length; i++)
          for (idim=0; idim<Ndim; idim++) 
            sten->Offset[i][idim] *= zone_coarseness;
      }
     
      if (Lcut_jac /*&& izone>0*/ && isten==U_ATTRACT){
         shorten_stencil(sten);
      }

      if (isten != POLYMER_CR){   /* take care of this later */
        renormalize_stencil(sten, vol_sten);
      }

      if (Iwrite == VERBOSE && Proc==0){
          print_out_stencil(isten, izone,icomp, jcomp, ifp);
      }

    } /* End of loop over components j */
    } /* End of loop over components i */
  } /* End of loop over stencil types */
  } /* End of loop over quadrature zones */


  /* add 2 nodes in each direction to compensate for inability
     to check beyond box boundaries when setting up Nodes_2_boundary array*/
  for (idim=0; idim<Ndim; idim++) Max_sten_length[idim] += 4;

  safe_free((void *) &el_offsets);
  safe_free((void *) &el_weights);

  if (Proc == 0) {
    fprintf(ifp2,"\n!!!!!!!!!!!! done setting up stencils !!!!!!!!!!!!!!!!!\n");
    fclose(ifp2);   

    if (Iwrite == VERBOSE) fclose(ifp);

    if (Iwrite != NO_SCREEN){
             printf("stencil setup took %g secs\n", MPI_Wtime()-t1);
             printf("---------------------------------------------------------------------------------\n");
    }
  }
}
/****************************************************************************/
void shorten_stencil(struct Stencil_Struct *sten)
{
/* eliminate small tail contributions from the
   U_ATTRACT stencil for Jacobian terms */

int i,j,k,length_tmp,count; 
double max_len;
struct Stencil_Struct sten_tmp;

length_tmp = sten->Length;

(&sten_tmp)->Offset = (int **) array_alloc(2, sten->Length, Ndim, sizeof(int));
(&sten_tmp)->Weight = (double *) array_alloc(1, sten->Length, sizeof(double));
if (Lhard_surf) (&sten_tmp)->HW_Weight = (double **)
      array_alloc(2, sten->Length, Nnodes_per_el_V, sizeof(double));

max_len=0.0;
for (i=0; i<sten->Length; i++){
   if (fabs(sten->Weight[i]) > max_len) max_len=fabs(sten->Weight[i]);
}

count=0;
for (i=0; i<sten->Length; i++){
   (&sten_tmp)->Weight[i] = sten->Weight[i];
   for (j=0; j<Ndim; j++)  (&sten_tmp)->Offset[i][j] = sten->Offset[i][j];
   if (Lhard_surf){
     for (j=0; j<Nnodes_per_el_V; j++)  
         (&sten_tmp)->HW_Weight[i][j] = sten->HW_Weight[i][j];
   }

   if (fabs(sten->Weight[i]) > max_len/Jac_threshold) count++;
}

safe_free((void *) &sten->Offset);
safe_free((void *) &sten->Weight);
if (Lhard_surf) safe_free((void *) &sten->HW_Weight);

sten->Offset = (int **) array_alloc(2, count, Ndim, sizeof(int));
sten->Weight = (double *) array_alloc(1, count, sizeof(double));
if (Lhard_surf) sten->HW_Weight = (double **)
      array_alloc(2, count, Nnodes_per_el_V, sizeof(double));

sten->Length=0;
for (i=0; i<length_tmp; i++){
   if (fabs((&sten_tmp)->Weight[i]) > max_len/Jac_threshold){
      sten->Weight[sten->Length] = (&sten_tmp)->Weight[i]; 
      for (j=0; j<Ndim; j++)
         sten->Offset[sten->Length][j] = (&sten_tmp)->Offset[i][j];

      if (Lhard_surf) {
           for (k=0; k<Nnodes_per_el_V; k++) 
                sten->HW_Weight[sten->Length][k] = (&sten_tmp)->HW_Weight[i][k];
      }
      sten->Length++;
   }
}

if (sten->Length != count){
   printf("problems with shortening the stencil\n");
   printf("number of entries: %d  expected number %d",sten->Length,count);
   exit(-1);
}

safe_free((void *) &(&sten_tmp)->Offset);
safe_free((void *) &(&sten_tmp)->Weight);
if (Lhard_surf) safe_free((void *) &(&sten_tmp)->HW_Weight);

return;
}
/****************************************************************************/

static void merge_load_stencil(struct Stencil_Struct *sten,
			       int **el_offsets, double *el_weights,
			       int *el_in_radius,int *index_sten)
/*
 * The variable sten in this routine is known globally as 
 * Stencil[isten][izone][icomp+Ncomp*jcomp]
 */
{
  int isten, j, k;
  int *el_offset_ptr;

  /* loop over element stencils stored in el_* arrays */
 
  for (j=0; j<Nnodes_per_el_V; j++) {

    el_offset_ptr = el_offsets[j];

    /* Move on to next entry in stencil if weight is zero */

    if (el_weights[j] != 0.0) {

       isten = ijk_to_isten_index(el_offset_ptr,el_in_radius);

        /* If the node is new, add to the end of the merge array */
       if (index_sten[isten] == -1) {
 
            sten->Weight[sten->Length] = el_weights[j];
            for (k=0; k<Ndim; k++) 
                 sten->Offset[sten->Length][k] = el_offset_ptr[k];
 
         /* If hard walls, init weights then set one for this local elem */
         if (Lhard_surf) {
           for (k=0; k<Nnodes_per_el_V; k++) 
                sten->HW_Weight[sten->Length][k] = 0.0;
           sten->HW_Weight[sten->Length][j] = el_weights[j];
         }
         index_sten[isten] = sten->Length;

         sten->Length++;
       }

       /* else the node is old, add the weights together */

       else {
         sten->Weight[index_sten[isten]] += el_weights[j];
         if (Lhard_surf)
           sten->HW_Weight[index_sten[isten]][j] = el_weights[j];
       }

    } /* End of if statement checking for non-zero weight */
  } /* End loop over all nodes in el_* arrays */

}
/****************************************************************************/
/* ijk_to_isten_index:  Given an i,j,k position of a node, return its node number */

int ijk_to_isten_index(int *ijk,int *el_in_radius)
{
  int ijk_tmp[3],nodes_plane=0,nodes_x[3],i,Nodes_x_tmp[3];
  for (i=0; i<Ndim; i++){
      Nodes_x_tmp[i] = round_to_int(Size_x[i]/Esize_x[i] + 1.)-1;
      if (Type_bc[i][0]==PERIODIC && (2*el_in_radius[i]+1) > Nodes_x_tmp[i]){
        nodes_x[i] = Nodes_x_tmp[i];
        if (ijk[i] < -Nodes_x_tmp[i]/2)
           ijk_tmp[i] = ijk[i] + Nodes_x_tmp[i];
        else if ( ijk[i] >= (Nodes_x_tmp[i]/2 + Nodes_x_tmp[i]%2) ) 
           ijk_tmp[i] = ijk[i] - Nodes_x_tmp[i];
        else ijk_tmp[i] = ijk[i];

        ijk_tmp[i] += Nodes_x_tmp[i]/2;
      }
      else{           
        nodes_x[i] = 2 * el_in_radius[i] + 1;
        ijk_tmp[i] = ijk[i] + el_in_radius[i];
      }
  }
  if (Ndim==3) nodes_plane = nodes_x[0]*nodes_x[1];

  switch (Ndim){
    case 1:   return (*ijk_tmp);
    case 2:   return (*ijk_tmp + ijk_tmp[1] * (*nodes_x));
    case 3:   return (*ijk_tmp + ijk_tmp[1] * (*nodes_x) + ijk_tmp[2] * nodes_plane);
  }
  return 0;
}
/****************************************************************************/

static void renormalize_stencil(struct Stencil_Struct *sten, double vol_sten)
{
   double sum = 0, ratio;
   int i,j;

   for (i=0; i < sten->Length; i++) 
     sum += sten->Weight[i];

   if (sum == vol_sten) return;

   if (fabs(sum - vol_sten) / vol_sten < 0.1) {
/*      if (Proc==0) printf("\tRenormalizing stencil from %g to known volume of %g\n",
            sum, vol_sten);*/
   }
   else if(Ndim !=3) {
      if(Proc==0) printf("WARNING: Stencil does not add up near correct volume %g %g\n",
             sum, vol_sten);
   }

   /* renormalize stencil */

   ratio = vol_sten/sum;
 
   for (i=0; i < sten->Length; i++) sten->Weight[i] *= ratio;

   if (Lhard_surf) {
     for (i=0; i < sten->Length; i++) 
       for (j=0; j < Nnodes_per_el_V; j++) 
         sten->HW_Weight[i][j] *= ratio;
   }
}
/****************************************************************************/

static void print_out_stencil(int isten, int izone,
			      int icomp, int jcomp, FILE *ifp)
{
  int i,j;
  struct Stencil_Struct *sten;
  double sum;

  /* set temporary stencil pointer, to avoid repeated indexing */

  sum=0.0;
  sten = &(Stencil[isten][izone][icomp + Ncomp*jcomp]);

  fprintf(ifp,"\nWriting out stencil for isten %d, izone %d, icomp %d, jcomp %d, Length %d\n\n",
              isten, izone, icomp, jcomp, sten->Length);

  fprintf(ifp," #  Node Offsets   Weight\n");
  for (i=0; i<sten->Length; i++) {
    if (Ndim<=3  && (sten->Length < 2500 || sten->Offset[i][0]==0)){
       fprintf(ifp,"%3d. ",i);
       for (j=0; j<Ndim; j++){ 
         fprintf(ifp,"  %3d  %7.4f",sten->Offset[i][j],Esize_x[j]*sten->Offset[i][j]);
       }
       fprintf(ifp,"  %e  ",sten->Weight[i]);
    }
    sum += sten->Weight[i];
    if (Lhard_surf &&(Ndim<=3 && (sten->Length < 2500 || sten->Offset[i][0]==0))){
         for (j=0; j<Nnodes_per_el_V; j++){
            fprintf(ifp," %5f",sten->HW_Weight[i][j]);
          }
    }
    if ( Ndim<3 || Ndim==3 &&sten->Length < 2500 || sten->Offset[i][0]==0) fprintf(ifp,"\n");
  }
  fprintf(ifp,"\tSUM OF WEIGHTS = %e\n\n",sum);

  fprintf(ifp,"============================================\n");

  return;
}

/****************************************************************************/

void sort2_int_double_array(int n, int ra[], double rb[],
			    int *rc[], int ndim)
/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array (JS)
*       modified to include int array and to number from 0
*
*       Sorts the array ra[1,..,n] in ascending numerical order using heapsort
*       algorityhm, while making the corresponding rearrangement of the
*       array rb[1,..,n] and rc[1,..,n][ndim].
*
*/

{
  int l,j,ir,i,idim;
  int rra, rrc[NDIM_MAX];
  double rrb;

  /*
   *  No need to sort if one or fewer items.
   */
  if (n <= 1) return;

  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra=ra[--l];
      rrb=rb[l];
      for (idim=0;idim<ndim;idim++) rrc[idim]=rc[l][idim];
    } else {
      rra=ra[ir];
      rrb=rb[ir];
      for (idim=0;idim<ndim;idim++) rrc[idim]=rc[ir][idim];
      ra[ir]=ra[1];
      rb[ir]=rb[1];
      for (idim=0;idim<ndim;idim++) rc[ir][idim]=rc[1][idim];
      if (--ir == 1) {
        ra[1]=rra;
        rb[1]=rrb;
        for (idim=0;idim<ndim;idim++) rc[1][idim]=rrc[idim];
        return;
      }
    }
    i=l;
    j=l << 1;
    while (j <= ir)     {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        rb[i]=rb[j];
        for (idim=0;idim<ndim;idim++) rc[i][idim]=rc[j][idim];
        j += (i=j);
      }
      else j=ir+1;
    }
    ra[i]=rra;
    rb[i]=rrb;
    for (idim=0;idim<ndim;idim++)  rc[i][idim]=rrc[idim];
  }
}

/****************************************************************************/

static void sort_stencil (struct Stencil_Struct *sten)
{
#define MAX_NODES_DIM 10000
  int i, *value;

  /* Assemble value array, with nodal offsets from current node */

  value = (int *) array_alloc(1, sten->Length, sizeof(int));

  for (i=0; i<sten->Length; i++) {
    value[i] = sten->Offset[i][0];
    if (Ndim>1) value[i] += sten->Offset[i][1] * MAX_NODES_DIM;
    if (Ndim>2) value[i] += sten->Offset[i][2] * MAX_NODES_DIM * MAX_NODES_DIM;
  }
  
  /* Sort Weight and Offset by value */
  /* Subtract 1 from array pointers since sort expects fortran indexing */
  sort2_int_double_array(sten->Length, value-1, sten->Weight-1,
                                          sten->Offset-1, Ndim);

  safe_free((void *) &value);
}

/****************************************************************************/
static double calc_sten_rad(int isten, int icomp, int jcomp)
{
  switch (isten) {
    case DELTA_FN:
        if (Sten_Type[POLYMER_CR]) return (Bond_ff[icomp][jcomp]);
        else                       return (HS_diam[icomp]/2.0);
    case THETA_FN:
        return (HS_diam[icomp]/2.0);
    case U_ATTRACT:
        return (Cut_ff[icomp][jcomp]);
    case THETA_CHARGE:
        return (Sigma_ff[icomp][icomp]);
    case POLYMER_CR:
        return (Cr_rad[icomp][jcomp]);
    case POLYMER_GAUSS:
        return (1.5*Sigma_ff[icomp][icomp]);   /* fix this later */
    case THETA_FN_SIG:
        return (Sigma_ff[icomp][icomp]);
    case DELTA_FN_BOND:
        return (Bond_ff[icomp][jcomp]);
  }
  return 0;
}
/****************************************************************************/
static double calc_sten_vol(int isten, int i, int j)
{
  double r_min,r_cut,vol_sten,r_max;
  switch (isten) {
    case DELTA_FN:
        if (Sten_Type[POLYMER_CR]) return (1.0);
        else                       return (PI * POW_DOUBLE_INT(HS_diam[i],2));
    case THETA_FN:
         return (PI * POW_DOUBLE_INT(HS_diam[i],3)/6.0);
    case U_ATTRACT:
        r_min = Sigma_ff[i][j] * pow(2.0,1.0/6.0);
        r_cut = Cut_ff[i][j];

        vol_sten =  (4.0/3.0)*PI*pow(r_min,3.0)*uLJatt_n_noshift(r_min,i,j)
                      - (4.0/3.0)*PI*pow(r_cut,3.0)*uLJatt_n_noshift(r_cut,i,j)
                      + uLJatt_n_int(r_cut,i,j) - uLJatt_n_int(r_min,i,j);
 
        return(vol_sten);

    case THETA_CHARGE:
        r_max = Sigma_ff[i][j];
        vol_sten = deltaC_MSA_int(r_max,i,j);

        return(vol_sten);

    case POLYMER_CR:  /* actually this is an integration of c(r) over all
                         volume. Maybe normalize this to the compressibility.
                         do this later  */
        vol_sten = 0.;
        return(vol_sten);

    case POLYMER_GAUSS:
        vol_sten = 1.;
        return(vol_sten);
  
    case THETA_FN_SIG:
/*        return (4.0 * PI * POW_DOUBLE_INT(Sigma_ff[i][i],3)/3.0);*/
/*        in order to avoid having to carry prefactors, we will set the 
          stencil volume to 1.0.  This will effectively perform the multiplication
          of 3/(4PI Sigma^3) times the native stencil.*/
          return (1.0);
          

    case DELTA_FN_BOND:
/*        return (4.0 * PI * POW_DOUBLE_INT(Bond_ff[i][j],2) );*/
/*        in order to avoid having to carry prefactors, we will set the 
          stencil volume to 1.0.  This will effectively perform the multiplication
          of 1/(4PI Sigma[icomp][jcomp]^2) times the native stencil.*/
          return (1.0);
  }
  return 0;
}
/****************************************************************************/
void set_gauss_quad(int ngp, double *gp, double *gw)
{

  if (ngp == 40) {
    gp[ 0] = ( 0.038772417506050821933 + 1.0)/2.0;  gw[ 0] = 0.077505947978424811264/2.0;
    gp[ 1] = ( 0.116084070675255208483 + 1.0)/2.0;  gw[ 1] = 0.077039818164247965588/2.0;
    gp[ 2] = ( 0.192697580701371099716 + 1.0)/2.0;  gw[ 2] = 0.076110361900626242372/2.0;
    gp[ 3] = ( 0.268152185007253681141 + 1.0)/2.0;  gw[ 3] = 0.074723169057968264200/2.0;
    gp[ 4] = ( 0.341994090825758473007 + 1.0)/2.0;  gw[ 4] = 0.072886582395804059061/2.0;
    gp[ 5] = ( 0.413779204371605001525 + 1.0)/2.0;  gw[ 5] = 0.070611647391286779695/2.0;
    gp[ 6] = ( 0.483075801686178712909 + 1.0)/2.0;  gw[ 6] = 0.067912045815233903826/2.0;
    gp[ 7] = ( 0.549467125095128202076 + 1.0)/2.0;  gw[ 7] = 0.064804013456601038075/2.0;
    gp[ 8] = ( 0.612553889667980237953 + 1.0)/2.0;  gw[ 8] = 0.061306242492928939167/2.0;
    gp[ 9] = ( 0.671956684614179548379 + 1.0)/2.0;  gw[ 9] = 0.057439769099391551367/2.0;
    gp[10] = ( 0.727318255189927103281 + 1.0)/2.0;  gw[10] = 0.053227846983936824355/2.0;
    gp[11] = ( 0.778305651426519387695 + 1.0)/2.0;  gw[11] = 0.048695807635072232061/2.0;
    gp[12] = ( 0.824612230833311663196 + 1.0)/2.0;  gw[12] = 0.043870908185673271992/2.0;
    gp[13] = ( 0.865959503212259503821 + 1.0)/2.0;  gw[13] = 0.038782167974472017640/2.0;
    gp[14] = ( 0.902098806968874296728 + 1.0)/2.0;  gw[14] = 0.033460195282547847393/2.0;
    gp[15] = ( 0.932812808278676533361 + 1.0)/2.0;  gw[15] = 0.027937006980023401098/2.0;
    gp[16] = ( 0.957916819213791655805 + 1.0)/2.0;  gw[16] = 0.022245849194166957262/2.0;
    gp[17] = ( 0.977259949983774262663 + 1.0)/2.0;  gw[17] = 0.016421058381907888713/2.0;
    gp[18] = ( 0.990726238699457006453 + 1.0)/2.0;  gw[18] = 0.010498284531152813615/2.0;
    gp[19] = ( 0.998237709710559200350 + 1.0)/2.0;  gw[19] = 0.004521277098533191258/2.0;
    gp[20] = (-0.038772417506050821933 + 1.0)/2.0;  gw[20] = gw[0];
    gp[21] = (-0.116084070675255208483 + 1.0)/2.0;  gw[21] = gw[1];
    gp[22] = (-0.192697580701371099716 + 1.0)/2.0;  gw[22] = gw[2];
    gp[23] = (-0.268152185007253681141 + 1.0)/2.0;  gw[23] = gw[3];
    gp[24] = (-0.341994090825758473007 + 1.0)/2.0;  gw[24] = gw[4];
    gp[25] = (-0.413779204371605001525 + 1.0)/2.0;  gw[25] = gw[5];
    gp[26] = (-0.483075801686178712909 + 1.0)/2.0;  gw[26] = gw[6];
    gp[27] = (-0.549467125095128202076 + 1.0)/2.0;  gw[27] = gw[7];
    gp[28] = (-0.612553889667980237953 + 1.0)/2.0;  gw[28] = gw[8];
    gp[29] = (-0.671956684614179548379 + 1.0)/2.0;  gw[29] = gw[9];
    gp[30] = (-0.727318255189927103281 + 1.0)/2.0;  gw[30] = gw[10];
    gp[31] = (-0.778305651426519387695 + 1.0)/2.0;  gw[31] = gw[11];
    gp[32] = (-0.824612230833311663196 + 1.0)/2.0;  gw[32] = gw[12];
    gp[33] = (-0.865959503212259503821 + 1.0)/2.0;  gw[33] = gw[13];
    gp[34] = (-0.902098806968874296728 + 1.0)/2.0;  gw[34] = gw[14];
    gp[35] = (-0.932812808278676533361 + 1.0)/2.0;  gw[35] = gw[15];
    gp[36] = (-0.957916819213791655805 + 1.0)/2.0;  gw[36] = gw[16];
    gp[37] = (-0.977259949983774262663 + 1.0)/2.0;  gw[37] = gw[17];
    gp[38] = (-0.990726238699457006453 + 1.0)/2.0;  gw[38] = gw[18];
    gp[39] = (-0.998237709710559200350 + 1.0)/2.0;  gw[39] = gw[19];
  }
  else if (ngp == 20) {
    gp[ 0] = ( 0.076526521133497333755 +1.0 )/2.0;  gw[ 0] = 0.152753387130725850698/2.0;
    gp[ 1] = ( 0.227785851141645078080 +1.0 )/2.0;  gw[ 1] = 0.149172986472603746788/2.0;
    gp[ 2] = ( 0.373706088715419560673 + 1.0)/2.0;  gw[ 2] = 0.142096109318382051329/2.0;
    gp[ 3] = ( 0.510867001950827098004 + 1.0)/2.0;  gw[ 3] = 0.131688638449176626898/2.0;
    gp[ 4] = ( 0.636053680726515025453 + 1.0)/2.0;  gw[ 4] = 0.118194531961518417312/2.0;
    gp[ 5] = ( 0.746331906460150792614 + 1.0)/2.0;  gw[ 5] = 0.101930119817240435037/2.0;
    gp[ 6] = (0.839116971822218823395 + 1.0)/2.0;  gw[ 6] = 0.083276741576704748725/2.0;
    gp[ 7] = (0.912234428251325905868 + 1.0)/2.0;  gw[ 7] = 0.062672048334109063570/2.0;
    gp[ 8] = (0.963971927277913791268 + 1.0)/2.0;  gw[ 8] = 0.040601429800386941331/2.0;
    gp[ 9] = (0.993128599185094924786 + 1.0)/2.0;  gw[ 9] = 0.017614007139152118312/2.0;
    gp[10] = (-0.076526521133497333755 +1.0 )/2.0;  gw[10] = gw[0];
    gp[11] = (-0.227785851141645078080 +1.0 )/2.0;  gw[11] = gw[1];
    gp[12] = (-0.373706088715419560673 + 1.0)/2.0;  gw[12] = gw[2];
    gp[13] = (-0.510867001950827098004 + 1.0)/2.0;  gw[13] = gw[3];
    gp[14] = (-0.636053680726515025453 + 1.0)/2.0;  gw[14] = gw[4];
    gp[15] = (-0.746331906460150792614 + 1.0)/2.0;  gw[15] = gw[5];
    gp[16] = (-0.839116971822218823395 + 1.0)/2.0;  gw[16] = gw[6];
    gp[17] = (-0.912234428251325905868 + 1.0)/2.0;  gw[17] = gw[7];
    gp[18] = (-0.963971927277913791268 + 1.0)/2.0;  gw[18] = gw[8];
    gp[19] = (-0.993128599185094924786 + 1.0)/2.0;  gw[19] = gw[9];
  }

  else if (ngp == 12) {
    gp[ 0] = ( 0.125233408511469 + 1.0)/2.0;  gw[ 0] = 0.249147045813403/2.0;
    gp[ 1] = ( 0.367831498998180 + 1.0)/2.0;  gw[ 1] = 0.233492536538355/2.0;
    gp[ 2] = ( 0.587317954286617 + 1.0)/2.0;  gw[ 2] = 0.203167426723066/2.0;
    gp[ 3] = ( 0.769902674194305 + 1.0)/2.0;  gw[ 3] = 0.160078328543346/2.0;
    gp[ 4] = ( 0.904117256370475 + 1.0)/2.0;  gw[ 4] = 0.106939325995318/2.0;
    gp[ 5] = ( 0.981560634246719 + 1.0)/2.0;  gw[ 5] = 0.047175336386512/2.0;
    gp[ 6] = (-0.125233408511469 + 1.0)/2.0;  gw[ 6] = 0.249147045813403/2.0;
    gp[ 7] = (-0.367831498998180 + 1.0)/2.0;  gw[ 7] = 0.233492536538355/2.0;
    gp[ 8] = (-0.587317954286617 + 1.0)/2.0;  gw[ 8] = 0.203167426723066/2.0;
    gp[ 9] = (-0.769902674194305 + 1.0)/2.0;  gw[ 9] = 0.160078328543346/2.0;
    gp[10] = (-0.904117256370475 + 1.0)/2.0;  gw[10] = 0.106939325995318/2.0;
    gp[11] = (-0.981560634246719 + 1.0)/2.0;  gw[11] = 0.047175336386512/2.0;
  }
  else if (ngp == 3) {
    gp[ 0] = ( 0.7745966692      + 1.0)/2.0;  gw[ 0] = 0.555555555555556/2.0;
    gp[ 1] = ( 0.0000000         + 1.0)/2.0;  gw[ 1] = 0.888888888888888/2.0;
    gp[ 2] = (-0.7745966692      + 1.0)/2.0;  gw[ 2] = 0.555555555555556/2.0;
  }
  else if (ngp == 1) {
    gp[ 0] = ( 0.0000000         + 1.0)/2.0;  gw[ 0] = 2.000000000000000/2.0;
  }
  else if (ngp == 5) {
    gp[ 0] = ( 0.906179845936884 + 1.0)/2.0;  gw[ 0] = 0.236926885056189/2.0;
    gp[ 1] = ( 0.538469310105683 + 1.0)/2.0;  gw[ 1] = 0.478628670499366/2.0;
    gp[ 2] = ( 0.00000           + 1.0)/2.0;  gw[ 2] = 0.568888888888890/2.0;
    gp[ 3] = (-0.538469310105683 + 1.0)/2.0;  gw[ 3] = 0.478628670499366/2.0;
    gp[ 4] = (-0.906179845936884 + 1.0)/2.0;  gw[ 4] = 0.236926885056189/2.0;
  }
  else if (ngp == 6) {
    gp[ 0] = ( 0.238619186083197 + 1.0)/2.0;  gw[ 0] = 0.467913934572691/2.0;
    gp[ 1] = ( 0.661209386466265 + 1.0)/2.0;  gw[ 1] = 0.360761573048139/2.0;
    gp[ 2] = ( 0.932469514203152 + 1.0)/2.0;  gw[ 2] = 0.171324492379170/2.0;
    gp[ 3] = (-0.238619186083197 + 1.0)/2.0;  gw[ 3] = 0.467913934572691/2.0;
    gp[ 4] = (-0.661209386466265 + 1.0)/2.0;  gw[ 4] = 0.360761573048139/2.0;
    gp[ 5] = (-0.932469514203152 + 1.0)/2.0;  gw[ 5] = 0.171324492379170/2.0;
  }
  else if (ngp == 9) {
    gp[ 0] = ( 0.324253423403809 + 1.0)/2.0;  gw[ 0] = 0.312347077040004/2.0;
    gp[ 1] = ( 0.613371432700590 + 1.0)/2.0;  gw[ 1] = 0.260610696402935/2.0;
    gp[ 2] = ( 0.836031107326636 + 1.0)/2.0;  gw[ 2] = 0.180648160694857/2.0;
    gp[ 3] = ( 0.968160239507626 + 1.0)/2.0;  gw[ 3] = 0.081274388361574/2.0;
    gp[ 4] = ( 0.000000000000000 + 1.0)/2.0;  gw[ 4] = 0.330239355001260/2.0;
    gp[ 5] = (-0.324253423403809 + 1.0)/2.0;  gw[ 5] = 0.312347077040004/2.0;
    gp[ 6] = (-0.613371432700590 + 1.0)/2.0;  gw[ 6] = 0.260610696402935/2.0;
    gp[ 7] = (-0.836031107326636 + 1.0)/2.0;  gw[ 7] = 0.180648160694857/2.0;
    gp[ 8] = (-0.968160239507626 + 1.0)/2.0;  gw[ 8] = 0.081274388361574/2.0;
  }
  else {
    printf("Requested Number of Gauss points not allowed (%d)\n",ngp);
    exit(-1);
  }

}
/****************************************************************************/

static double get_weight_from_stencil(int isten, int icomp, int jcomp, double rsq,
				      double R, int ngpu, double *gpu, double *gwu)
{
  double temp, zmax, z, rho, rmin, rlast_nz, r_upp,slope_dr,r_low,zsq,rx_low;
  int i, irmin;
  double sigsq;

  switch (isten) {

    case DELTA_FN:
    case DELTA_FN_BOND:
          if (!Sten_Type[POLYMER_CR]){
            if (Ndim == 1)       return(2.0 * PI * R);
            else if (Ndim == 2)  return(2.0 / sqrt(1.0-rsq));
            else                 return(1.0);
          }
	  /* ALF: fixed normalization of delta funcs */
          else{
	    sigsq = Sigma_ff[icomp][jcomp]*Sigma_ff[icomp][jcomp];
            if (Ndim == 1)       return(R / (2.0*sigsq));
            else if (Ndim == 2)  return(0.5 / (sqrt(1.0-rsq) * PI*sigsq));
            else                 return(0.25 / (PI*sigsq));
          }

    case THETA_FN:
    case THETA_FN_SIG:
            if (Ndim == 1)       return((1.0-rsq) * PI*R*R);
            else if (Ndim == 2)  return(sqrt(1.0-rsq) * 2.0*R);
            else                 return(1.0);

    case U_ATTRACT:

            if (Ndim == 1) {
               temp = 0.0;
               zmax = sqrt(1 - rsq);
               for (i=0; i < ngpu; i++) {

                  z = zmax * gpu[i];
                  rho = sqrt(rsq + z*z) * Cut_ff[icomp][jcomp];
                  temp += gwu[i] * z * uLJatt_n(rho, icomp, jcomp);
                  if (uLJatt_n(rho,icomp,jcomp) > 0.)
                       printf("rij: %f  uLJ: %f\n",rho,
                              uLJatt_n(rho,icomp,jcomp));
               }
               return(2.0 * PI * temp * Cut_ff[icomp][jcomp]
                                      * Cut_ff[icomp][jcomp] * zmax);
            }
            else if (Ndim == 2) {
               temp = 0.0;
               zmax = sqrt(1 - rsq);
               for (i=0; i < ngpu; i++) {
                  z = zmax * gpu[i];
                  rho = sqrt(rsq + z*z) * Cut_ff[icomp][jcomp];
                  temp += gwu[i] * uLJatt_n(rho, icomp, jcomp);
               }
               return(2.0 * temp * Cut_ff[icomp][jcomp] * zmax);
            }
            else {
              rho = sqrt(rsq) * Cut_ff[icomp][jcomp];
              temp = uLJatt_n(rho, icomp, jcomp);
              return(temp);
            }

    case THETA_CHARGE:
            if (Ndim == 1) {
               temp = 0.0;
               zmax = sqrt(1 - rsq);
               for (i=0; i < ngpu; i++) {
                  z = zmax * gpu[i];
                  /* perhaps rho doesn't need Sigma_ff multiplier */
                  rho = sqrt(rsq + z*z) * Sigma_ff[icomp][jcomp];
                  temp += gwu[i] * z *  deltaC_MSA(rho, icomp, jcomp);
               }
               return(2.0 * PI * temp * Sigma_ff[icomp][jcomp]
                                      * Sigma_ff[icomp][jcomp] * zmax);
            }
            else if (Ndim == 2) {
               temp = 0.0;
               zmax = sqrt(1 - rsq);
               for (i=0; i < ngpu; i++) {
                  z = zmax * gpu[i];
                  /* perhaps rho doesn't need Sigma_ff multiplier */
                  rho = sqrt(rsq + z*z) * Sigma_ff[icomp][jcomp];
                  temp += gwu[i] *  deltaC_MSA(rho, icomp, jcomp);
               }
               return(2.0 * temp * Sigma_ff[icomp][jcomp] * zmax);
            }
            else {
              rho = sqrt(rsq) * Sigma_ff[icomp][jcomp];
              temp = deltaC_MSA(rho, icomp, jcomp);
              return(temp);
            }

    case POLYMER_CR:  
        if (Ndim < 3) {
          rmin = rsq;
          zsq = rmin*rmin;
          irmin = (int) (rmin/Deltar_cr + 0.00000001);
     /*   last_nz_cr = (int) (R/Deltar_cr + 0.00000001);     */
          rlast_nz   = Deltar_cr*Last_nz_cr;          
          temp = 0.;
          if (irmin < Last_nz_cr) {
            r_upp = Deltar_cr*(irmin+1);
            if (r_upp > R) r_upp = R;
            slope_dr = (Rism_cr[icomp][jcomp][irmin+1] - 
                       Rism_cr[icomp][jcomp][irmin]);
            rx_low = 0.;
            temp = int_cr(rmin,r_upp,slope_dr,icomp,jcomp,irmin,zsq,&rx_low);
            for (irmin = irmin+1; irmin < Last_nz_cr; irmin++) {
             if (r_upp < R){
              r_low = r_upp;
              r_upp = Deltar_cr*(irmin+1);
              if (r_upp > R) r_upp = R;
              slope_dr = (Rism_cr[icomp][jcomp][irmin+1] - 
                       Rism_cr[icomp][jcomp][irmin]);
             temp += int_cr(r_low,r_upp,slope_dr,icomp,jcomp,irmin,zsq,&rx_low);
             }
            }
          }
          if (rlast_nz < R) {
            if (rlast_nz > rmin){
               r_low = rlast_nz;
               rx_low = sqrt(r_low*r_low - zsq);
            }
            else{
               r_low = rmin;
               rx_low = 0.;
            }
            r_upp = R;
            slope_dr = (Rism_cr[icomp][jcomp][Last_nz_cr] - 
                       Rism_cr[icomp][jcomp][Last_nz_cr-1]);
        temp += int_cr(r_low,r_upp,slope_dr,icomp,jcomp,Last_nz_cr,zsq,&rx_low);
          }          
          return(temp);
        }
        else {
          rmin = rsq;
          irmin = (int) (rmin/Deltar_cr + 0.00000001);
          if (irmin < Last_nz_cr) 
            slope_dr = (Rism_cr[icomp][jcomp][irmin+1] - 
                       Rism_cr[icomp][jcomp][irmin]);
          else
            slope_dr = (Rism_cr[icomp][jcomp][irmin] - 
                       Rism_cr[icomp][jcomp][irmin-1]);
          temp = Rism_cr[icomp][jcomp][irmin] + 
                                   slope_dr*(rmin - irmin*Deltar_cr)/Deltar_cr;
          return(temp);
        }

    case POLYMER_GAUSS:   /* fix the cut-off values here */

            r_upp = -1.5*R*R/(Gauss_a*Gauss_a);
            if (Ndim == 1) {
               temp = 0.0;
               zmax = sqrt(1 - rsq);
               for (i=0; i < ngpu; i++) {

                  z = zmax * gpu[i];
                  rho = (rsq + z*z);
                  temp += gwu[i] * z * exp(rho*r_upp);
/*                temp += gwu[i] * z * gauss(rho, icomp, jcomp);     */
               }
               return(2.0 * PI * temp * R*Gauss_k * R*zmax);
            }
            else if (Ndim == 2) {
               temp = 0.0;
               zmax = sqrt(1 - rsq);
               for (i=0; i < ngpu; i++) {
                  z = zmax * gpu[i];
                  rho = (rsq + z*z);
                  temp += gwu[i] * exp(rho*r_upp);
               }
               return(2.0 * temp * R*Gauss_k * zmax);
            }
            else {
              temp = exp(rsq*r_upp);
              return(temp*Gauss_k);
            }
  }
  return 0;
}

/*******************************************************************************
  int_cr:  integrate the direct correlation function                         */
double int_cr(double r_low,double r_upp,double slope_dr,int icomp,int jcomp,
                  int irmin, double zsq, double *rx_low)
{
  double temp,rusq,rlsq,rx_upp;
/*
if (icomp==0 && jcomp==0 && fabs(r_low-3.411270)<1.e-4 && fabs(r_upp- 3.42500)<1.e-4){
   printf("Deltar_cr=%9.6f\n",Deltar_cr);
}
*/
  rusq = r_upp*r_upp;
  rlsq = r_low*r_low;
  temp=0.0;
  if (Ndim == 1)   {
    temp = 2.*PI*(
      ( rusq*(Rism_cr[icomp][jcomp][irmin] - slope_dr*irmin)/2. +
        slope_dr*rusq*r_upp/(3.*Deltar_cr) ) -
      ( rlsq*(Rism_cr[icomp][jcomp][irmin] - slope_dr*irmin)/2. +
        slope_dr*rlsq*r_low/(3.*Deltar_cr) )   );
  }
  else{
    rx_upp = sqrt(r_upp*r_upp - zsq);
    temp = 2.*(rx_upp - *rx_low)*(Rism_cr[icomp][jcomp][irmin] - slope_dr*irmin)
           + slope_dr*( (r_upp*rx_upp - r_low* *rx_low) +
           zsq*log((r_upp + rx_upp)/(r_low + *rx_low)) )/Deltar_cr;
    *rx_low = rx_upp;
  }
  return temp;
}
/*******************************************************************************
  gauss:      calculate Gaussian function                         */
double gauss(double r, int i, int j)
{
  void gcf(double *gammcf, double a, double psq, double *gln);
  void gser(double *gamser, double a, double psq, double *gln);

  double psq,gamser,gammcf,gln,a;

  a = 0.5;
  psq = 1.5*(r/Gauss_a)*(r/Gauss_a);
  if (psq < 1.5) {
    gser(&gamser,a,psq,&gln);
    return gamser;
  }
  else {
    gcf(&gammcf,a,psq,&gln);
    return 1.0 - gammcf;
  }
}
        

#define itmax 100
#define eps 3.e-7

void gser(double *gamser, double a, double x, double *gln)
{
  double gammln(double xx);
  int n;
  double sum,del,ap;

  *gln = gammln(a);
    if (x <= 0.) {
        *gamser=0.;
        return;
    }
    else {
      ap=a;
      del=sum=1./a;
      for (n=1; n<=itmax; n++) {
        ++ap;
        del *= x/ap;
        sum += del;
        if (fabs(del) < fabs(sum)*eps) {
          *gamser=sum*exp(-x+a*log(x)-(*gln));
          return;
        }
      }
      printf("a too large, itmax too small in gcf\n");
      return;
    }
}

#define itmax 100
#define eps 3.e-7
#define fpmin 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
  double gammln(double xx);
  int i;
  double an,b,c,d,del,h;

  *gln = gammln(a);
  b = x + 1.0 - a;
  c = 1.0/fpmin;
  d = 1.0/b;
  h = d;
  for (i=1; i<itmax; i++) {
    an = -i*(i-a);
    b += 2.0;
    d = an*d + b;
    if (fabs(d) < fpmin) d = fpmin;
    c = b + an/c;
    if (fabs(c) < fpmin) c = fpmin;
    d = 1.0/d;
    del = d*c;
    h *= del;
    if (fabs(del - 1.0) < eps) break;
  }
  if (i > itmax) printf("a too large, itmax too small in gcf\n");
  *gammcf = exp(-x + a*log(x) - (*gln))*h;
}


double gammln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009173e0,-86.50532033e0,24.01409822e0,
          -1.231739516e0,.120858003e-2,-.536382e-5};
  int j;

      y = x = xx;
      tmp = x + 5.5;
      tmp -= (x+0.5)*log(tmp);
      ser= 1.000000000190015;
      for (j=0; j <= 5; j++) ser += cof[j]/++y;
      return -tmp+log(2.50662827465e0*ser/x);
}
  
/****************************************************************************/
static int calc_in_out_on(double *x_l, double *x_r, double sten_rad)
{
  int inflag = FALSE, outflag = FALSE;
  double r;

  switch (Ndim) {

    case 1:
       /* check 2 end nodes */

       r = fabs(x_l[0]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;
       r = fabs(x_r[0]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;

       break;

    case 2:
       /* check 4 corner nodes */

       r = sqrt(x_l[0]*x_l[0] + x_l[1]*x_l[1]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;
       r = sqrt(x_l[0]*x_l[0] + x_r[1]*x_r[1]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;
       r = sqrt(x_r[0]*x_r[0] + x_l[1]*x_l[1]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;
       r = sqrt(x_r[0]*x_r[0] + x_r[1]*x_r[1]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;

       break;

    case 3:
       /* check 8 corner nodes */

       r = sqrt(x_l[0]*x_l[0] + x_l[1]*x_l[1] + x_l[2]*x_l[2]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;
       r = sqrt(x_l[0]*x_l[0] + x_r[1]*x_r[1] + x_l[2]*x_l[2]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;
       r = sqrt(x_r[0]*x_r[0] + x_l[1]*x_l[1] + x_l[2]*x_l[2]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;
       r = sqrt(x_r[0]*x_r[0] + x_r[1]*x_r[1] + x_l[2]*x_l[2]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;
       r = sqrt(x_l[0]*x_l[0] + x_l[1]*x_l[1] + x_r[2]*x_r[2]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;
       r = sqrt(x_l[0]*x_l[0] + x_r[1]*x_r[1] + x_r[2]*x_r[2]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;
       r = sqrt(x_r[0]*x_r[0] + x_l[1]*x_l[1] + x_r[2]*x_r[2]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;
       r = sqrt(x_r[0]*x_r[0] + x_r[1]*x_r[1] + x_r[2]*x_r[2]);
         if (r > sten_rad) outflag = TRUE; else inflag = TRUE;

       break;
  }

  /* return 1 if the element is completely outside, -1 inside, 0 stradling */
  if      (inflag == FALSE && outflag == TRUE ) return(1);
  else if (inflag == TRUE  && outflag == FALSE) return(-1);
  else return(0);

}
/****************************************************************************/
