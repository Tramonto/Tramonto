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
 *  FILE: dft_stencil.c
 *
 *  This file contains the routines to precalculate the "stencils,"
 *  which are the integration schemes for the non-local interactions.
 *  The results will be strored in the global variable Stencil, which
 *  is a 2D array of Stencil_Struct (defined in dft_stencil_const.h).
 *
 */

#include "dft_stencil.h"

void calc_stencils(void)
{
 /*
  * Local variable declarations
  */

  char *yo = "calc_stencils";
  int izone, isten, istart,iend,jstart,jend,kstart,kend,
      icomp, i, j, k;   /* counters in for loops         */
  int jcomp, imax, jmax, idim, iii; /* counters in for loops      */

  int    **el_offsets;
  double *el_weights;

  int max_sten_nodes, el_in_radius[3], in_out_on_flag;
  double vol_sten, sten_rad,t1;
  double x_left[3], x_right[3],x_mid[3],r_center_sq;

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
    if( (ifp2 = fopen("dft_out.lis","a+")) ==NULL){
      printf("Can't open file dft_out.lis\n");
      exit(1);
    }
  }

 /*
  * Allocate the Stencil variable to be a 3D array of structures
  * Now it can be accessed as Stencil[i][j][k].Length
  */

  Stencil = (struct Stencil_Struct ***) array_alloc
             (3, NSTEN, Nzone, Ncomp*Ncomp, sizeof(struct Stencil_Struct));

  el_offsets = (int **) array_alloc (2, Nnodes_per_el_V, Ndim, sizeof(int));
  el_weights = (double *) array_alloc (1, Nnodes_per_el_V, sizeof(double));

  /* zero the maximum stencil length .... used later for
     setting up local coordinate systems */
  for (idim=0; idim<Ndim; idim++)  Max_sten_length[idim] = 0;
 
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

                           /* for case of general MSA - we want to precalculate certain parameters
                              outside of the gauss quadrature loops needed for stencil set up. */
      if (isten==THETA_CR_GENERAL_MSA){
         printf("isten=%d calling precalc_GENmsa_params\n",isten);
          precalc_GENmsa_params(Rho_b,X_MSA,N_MSA,Gamma_MSA);
      }

                           /* set index if this stencil is first order (depends on icomp only) or 
                              second order (depends on icomp,jcomp pair) */
    imax = Ncomp;
    jmax = stencil_Njcomp_switch(isten);

    for (icomp=0; icomp<imax; icomp++) {
    for (jcomp=0; jcomp<jmax; jcomp++) {

      /* calculate stencil volume (e.g. 4/3 pi r^3), radius, and element vol */

      vol_sten = stencil_volume_switch(isten, icomp, jcomp);
      sten_rad = stencil_radius_switch(isten, icomp, jcomp);

      /*
       * Figure out elements that form a bounding box for the stencil 
       * and the maximum number of nodes in the stencil
       */

      max_sten_nodes = 1.0;
      el_in_radius[0]=0; el_in_radius[1]=0; el_in_radius[2]=0;
      for (i=0; i<Ndim; i++) {
        el_in_radius[i] = (int) ((float)sten_rad/(float) esize_zone[i] + 0.999);
        max_sten_nodes *= 2 * el_in_radius[i] + 1;
        if (el_in_radius[i]*zone_coarseness > Max_sten_length[i]) 
              Max_sten_length[i] = el_in_radius[i]*zone_coarseness;
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
        
      /* Loop through all elements in bounding box, and set up stencil */

                     istart=-el_in_radius[0]; iend=el_in_radius[0];
      if (Ndim >=2){ jstart=-el_in_radius[1]; jend=el_in_radius[1]; }
      else         { jstart=0; jend=1;}
      if (Ndim==3) { kstart=-el_in_radius[2]; kend=el_in_radius[2]; }
      else         { kstart=0; kend=1;}


      for (i = istart; i < iend; i++) {
         for (j = jstart; j < jend; j++) {
            for (k = kstart; k < kend; k++) {

             setup_stencil_offsets(i,j,k,el_offsets);
             switch(Ndim){
                case 3:
                    x_left[2]  = k * esize_zone[2];
                    x_right[2] = (k+1) * esize_zone[2];
                case 2:
                     x_left[1]  = j * esize_zone[1];
                     x_right[1] = (j+1) * esize_zone[1];
                case 1:
                     x_left[0]  = i * esize_zone[0];
                     x_right[0] = (i+1) * esize_zone[0];
             }
             for (idim=0;idim<Ndim;idim++) x_mid[idim]=x_left[idim] + 0.5*esize_zone[idim];
             r_center_sq=0.0;
             for (idim=0;idim<Ndim;idim++) r_center_sq+=x_mid[idim]*x_mid[idim];

             for (iii=0; iii < Nnodes_per_el_V; iii++) el_weights[iii] = 0.0;

             /* test if element is out, in, or stradles radius */
             in_out_on_flag = calc_in_out_on(x_left, x_right, sten_rad);


             if ( in_out_on_flag == 1 || 
                   (in_out_on_flag == -1 && Ndim==3 && stencil_deltaLogical(isten) )){
                        /* nothing to do -- element is either outside of stencil radius or is completely
                           inside a stencil defined by a delta function in 3-dimensions. */ 
             } 
             else{  /* need to get some contributions for the stencil */
                 if (in_out_on_flag == -1) { /* element completely within stencil radius - do Gauss quadrature */


                     calc_sten_weight_el_in_radius(isten,sten_rad,icomp,jcomp,r_center_sq, 
                                                            x_left,esize_zone,el_weights);

                 }
                 else{   /* element stradles stencil boundary - do lots of dumb quadrature*/

                     calc_sten_weight_el_on_boundary(isten,sten_rad,icomp,jcomp,r_center_sq, 
                                                             x_left,esize_zone,el_weights);

                 }
                 merge_load_stencil(sten, el_offsets, el_weights,el_in_radius,index_sten);
             }
 

          } /* end loop in k(z) */
        }   /* end loop in j (y) */
      }  /* end loop in i (x) */

      safe_free((void *)&index_sten);

      /* coarse stencils for big esize_zone must be scaled to real Esize*/

      if (izone > 0){
        for (i=0; i<sten->Length; i++)
          for (idim=0; idim<Ndim; idim++) 
            sten->Offset[i][idim] *= zone_coarseness;
      }
    
       
      if (Lcut_jac /*&& izone>0*/ && isten==THETA_PAIRPOT_RCUT){
         shorten_stencil(sten);
      }

      if ( !(fabs(vol_sten-(double)NO_RENORMALIZATION_FLAG) < 1.e-8) ) {
            if (Iwrite==VERBOSE && Proc==0){
            printf("Renormalizing stencil: isten=%d",isten);
            switch(isten){
                case DELTA_FN_R: printf(" (DELTA_FN_R stencil)\n"); break;
                case THETA_FN_R: printf(" (THETA_FN_R stencil)\n"); break;
                case THETA_PAIRPOT_RCUT: printf(" (THETA_PAIRPOT_RCUT stencil)\n"); break;
                case THETA_CR_RPM_MSA: printf(" (THETA_CR_RPM_MSA stencil)\n"); break;
                case THETA_CR_DATA: printf(" (THETA_CR_DATA stencil)\n"); break;
                case THETA_FN_SIG: printf(" (THETA_FN_SIG stencil)\n"); break;
                case DELTA_FN_BOND: printf(" (DELTA_FN_BOND stencil)\n"); break;
                case THETA_CR_GENERAL_MSA: printf(" (THETA_CR_GENERAL_MSA stencil)\n"); break;
                default: printf("problem with isten switch exititng\n"); exit(-1); break;
            }
            }
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
  return;
}
/****************************************************************************/
void setup_stencil_offsets(int i, int j, int k, int **el_offsets)
{
   switch (Ndim) {

      case 3:
           el_offsets[4][0] = i;
           el_offsets[5][0] = i + 1;
           el_offsets[6][0] = i;
           el_offsets[7][0] = i + 1;

           el_offsets[4][1] = j;
           el_offsets[5][1] = j;
           el_offsets[6][1] = j + 1;
           el_offsets[7][1] = j + 1;

           el_offsets[0][2] = k;
           el_offsets[1][2] = k;
           el_offsets[2][2] = k;
           el_offsets[3][2] = k;
           el_offsets[4][2] = k + 1;
           el_offsets[5][2] = k + 1;
           el_offsets[6][2] = k + 1;
           el_offsets[7][2] = k + 1;

      case 2:
           el_offsets[2][0] = i;
           el_offsets[3][0] = i + 1;

           el_offsets[0][1] = j;
           el_offsets[1][1] = j;
           el_offsets[2][1] = j + 1;
           el_offsets[3][1] = j + 1;

      case 1:
           el_offsets[0][0] = i;
           el_offsets[1][0] = i + 1;
      }
      return;
}
/****************************************************************************/
void calc_sten_weight_el_in_radius(int isten, double sten_rad, int icomp,
                                   int jcomp, double r_center_sq, double *x_left,
                                   double *esize_zone,double *el_weights)
{
   int ig=0,jg=0,kg=0;
   double point[3],weight,inv_sten_rad_sq,el_vol,radius_sq;
   int ngp, ngpu,i;
   double gp[12],gw[12],gpu[40], gwu[40];

   inv_sten_rad_sq = 1.0/(sten_rad*sten_rad);
   el_vol=1.0; for (i=0; i<Ndim; i++) el_vol *= esize_zone[i];

                                    /* stencil quadrature */
   ngp=stencil_quadGauss_switch(isten,r_center_sq);
                                    /* integrand quadrature off-dimensions*/
   ngpu=stencil_quadGaussIntegrand_switch(isten,r_center_sq);

   set_gauss_quad(ngp, gp, gw);
   if (ngpu>0) set_gauss_quad(ngpu, gpu, gwu);

   for (ig=0; ig < ngp; ig++) {
      point[0] = x_left[0] + gp[ig] * esize_zone[0];
      switch(Ndim){
         case 3:
             for (jg=0; jg < ngp; jg++) {
                point[1] = x_left[1] + gp[jg] * esize_zone[1];

                for (kg=0; kg < ngp; kg++) {
                   point[2] = x_left[2] + gp[kg] * esize_zone[2];

                   radius_sq = (point[0]*point[0] + point[1]*point[1]
                                 + point[2]*point[2])*inv_sten_rad_sq;

                   weight = stencil_GetWeight_switch(isten, icomp, jcomp, 
                                        radius_sq,sten_rad, ngpu, gpu, gwu);

                   weight *= gw[ig] * gw [jg] * gw[kg] * el_vol;

                   distribute_weight_to_nodes(weight,gp,el_weights,ig,jg,kg);
 
                }
             }
             break;

         case 2:
             for (jg=0; jg < ngp; jg++) {
                   point[1] = x_left[1] + gp[jg] * esize_zone[1];

                   radius_sq = (point[0]*point[0] + point[1]*point[1]) *inv_sten_rad_sq;
                   weight = stencil_GetWeight_switch(isten, icomp, jcomp, 
                                        radius_sq,sten_rad, ngpu, gpu, gwu);

                   weight *= gw[ig] * gw [jg] * el_vol;
                   distribute_weight_to_nodes(weight,gp,el_weights,ig,jg,kg);
              
             }
             break;

         case 1:
             radius_sq = point[0]*point[0]*inv_sten_rad_sq;
             weight = stencil_GetWeight_switch(isten, icomp, jcomp, 
                                               radius_sq,sten_rad, ngpu, gpu, gwu);
             weight *= gw[ig] * el_vol;
             distribute_weight_to_nodes(weight,gp,el_weights,ig,jg,kg);

             break;
  
         }
    }
    return;
}
/*********************************************************************************************/
void calc_sten_weight_el_on_boundary(int isten,double sten_rad,int icomp,int jcomp,double r_center_sq, 
                                double *x_left,double *esize_zone,double *el_weights)
{
   int npt,ngpu;
   double inv_npt,el_vol,ok_distance=0.0,point[3],radius_sq,inv_sten_rad_sq,tmp_sum;
   double gpu[40], gwu[40],qp[3],weight;
   int ig=0,jg=0,kg=0,i,idim;

   inv_sten_rad_sq = 1.0/(sten_rad*sten_rad);
   el_vol=1.0; for (i=0; i<Ndim; i++) el_vol *= esize_zone[i];

                                    /* stencil quadrature */
   npt = stencil_quadBoundaryEl_switch(isten);
                                    /* integrand quadrature off-dimensions*/
   ngpu=stencil_quadGaussIntegrand_switch(isten,r_center_sq);
   if (ngpu>0) set_gauss_quad(ngpu, gpu, gwu);


   inv_npt = 1.0 / (double) npt;

   tmp_sum=0.0;
   if (Ndim==3 && stencil_deltaLogical(isten)){       
       for (idim=0;idim<Ndim;idim++) tmp_sum +=esize_zone[idim]*esize_zone[idim];
       ok_distance=0.5*inv_npt*sqrt(tmp_sum);
   }

   for (ig=0; ig < npt; ig++) {
      qp[0] =  ((double)ig + 0.5) * inv_npt;
      point[0] = x_left[0] + qp[0] * esize_zone[0];

      switch(Ndim){
        case 1:
           radius_sq = point[0]*point[0]*inv_sten_rad_sq;

           if (radius_sq <= 1.000001) { /* only use gauss points on the delta function boundary */
              weight = stencil_GetWeight_switch(isten, icomp, jcomp, 
                                       radius_sq,sten_rad, ngpu, gpu, gwu);
              if (fabs(radius_sq-1.0) <= 1.e-8) weight *= 0.5;
              weight *= inv_npt * el_vol; 

              distribute_weight_to_nodes(weight,qp,el_weights,0,1,2);
           }
           break;
  
        case 2:
           for (jg=0; jg < npt; jg++) {

               qp[1] =  ((double)jg + 0.5) * inv_npt;
               point[1] = x_left[1] + qp[1] * esize_zone[1];
               radius_sq = (point[0]*point[0] + point[1]*point[1])
                           * inv_sten_rad_sq;

               if (radius_sq <= 1.0000001) {
                  weight = stencil_GetWeight_switch(isten, icomp, jcomp, 
                                     radius_sq,sten_rad, ngpu, gpu, gwu);

                  if (fabs(radius_sq-1.0) <= 1.e-8) weight *= 0.5;
                  weight *= inv_npt * inv_npt * el_vol; 
                  distribute_weight_to_nodes(weight,qp,el_weights,0,1,2);
               }
            }
            break;

        case 3:
             for (jg=0; jg < npt; jg++) {

               qp[1] =  ((double)jg + 0.5) * inv_npt;
               point[1] = x_left[1] + qp[1] * esize_zone[1];

               for (kg=0; kg < npt; kg++) {

                 qp[2] =  ((double)kg + 0.5) * inv_npt;
                 point[2] = x_left[2] + qp[2] * esize_zone[2];
                 radius_sq = (point[0]*point[0] + point[1]*point[1]
                               + point[2]*point[2])*inv_sten_rad_sq;

                 if ( ((radius_sq <= 1.00000001) && !stencil_deltaLogical(isten)) || 
                      ((stencil_deltaLogical(isten)) && (fabs(radius_sq-1.0) < 2*ok_distance)) ) {

                     weight = stencil_GetWeight_switch(isten, icomp, jcomp, 
                                    radius_sq,sten_rad, ngpu, gpu, gwu);

                     if (fabs(radius_sq-1.0) <= 1.e-8 && !stencil_deltaLogical(isten)) weight *= 0.5;
                     weight *=  (inv_npt*inv_npt*inv_npt) * el_vol; 

                     distribute_weight_to_nodes(weight,qp,el_weights,0,1,2);
                 }
               }
             }
             break;
      }  /* end of Ndim switch */

   } /* end of ig loop*/
   return;
}

/****************************************************************************/
void distribute_weight_to_nodes(double weight,double *qp,double *el_weights,int ig, int jg, int kg)
{
   switch(Ndim){
      case 3:
       el_weights[0] += weight * (1.0-qp[ig]) * (1.0-qp[jg]) * (1.0-qp[kg]);
       el_weights[1] += weight * qp[ig]  * (1.0-qp[jg]) * (1.0-qp[kg]);
       el_weights[2] += weight * (1.0-qp[ig]) *      qp[jg]  * (1.0-qp[kg]);
       el_weights[3] += weight * qp[ig] *       qp[jg]  * (1.0-qp[kg]); 
       el_weights[4] += weight * (1.0-qp[ig]) * (1.0-qp[jg]) *      qp[kg];
       el_weights[5] += weight * qp[ig]  * (1.0-qp[jg]) *      qp[kg];
       el_weights[6] += weight * (1.0-qp[ig]) *      qp[jg]  *      qp[kg];
       el_weights[7] += weight * qp[ig]  *      qp[jg]  *      qp[kg]; 
       break;

      case 2:
       el_weights[0] += weight * (1.0-qp[ig]) * (1.0-qp[jg]); 
       el_weights[1] += weight * qp[ig] * (1.0-qp[jg]); 
       el_weights[2] += weight * (1.0-qp[ig]) * qp[jg]; 
       el_weights[3] += weight * qp[ig] * qp[jg]; 
       break;

      case 1:
        el_weights[0] += weight * (1.0-qp[ig]); 
        el_weights[1] += weight * qp[ig]; 
        break;
   }
   return;
}
/****************************************************************************/
void shorten_stencil(struct Stencil_Struct *sten)
{
/* eliminate small tail contributions from the
   THETA_PAIRPOT_RCUT stencil for Jacobian terms */

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

void merge_load_stencil(struct Stencil_Struct *sten,
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
           sten->HW_Weight[index_sten[isten]][j] += el_weights[j];
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

void renormalize_stencil(struct Stencil_Struct *sten, double vol_sten)
{
   double sum = 0.0, ratio;
   int i,j;

   for (i=0; i < sten->Length; i++) sum += sten->Weight[i];

   if (Iwrite==VERBOSE && Proc==0) printf("\t before normalization vol_sten=%g  sum_sten=%g\n",vol_sten,sum);

   if (sum == vol_sten) return;

   /* renormalize stencil */

   ratio = vol_sten/sum;
 
   for (i=0; i < sten->Length; i++) sten->Weight[i] *= ratio;

   if (Lhard_surf) {
     for (i=0; i < sten->Length; i++) 
       for (j=0; j < Nnodes_per_el_V; j++) 
         sten->HW_Weight[i][j] *= ratio;
   }

   sum=0.0;
   for (i=0; i < sten->Length; i++) sum += sten->Weight[i];
   if (Iwrite==VERBOSE && Proc==0) printf("\t after normalization vol_sten=%g  sum_sten=%g\n",vol_sten,sum);
}
/****************************************************************************/

void print_out_stencil(int isten, int izone,
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
int calc_in_out_on(double *x_l, double *x_r, double sten_rad)
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
