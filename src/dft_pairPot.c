/*
//@HEADER
// ********************************************************************
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/*
 *  FILE: dft_pairPot.c
 *
 *  This file contains routines general to any pair potential.
 *
 */

#include "dft_pairPot.h"


/*******************************************************************************************/
/*setup_pairPot:  some logic needed to set up pair potentials including locting potential minima and zeros */
void setup_pairPotentials(char *file_echoinput){
  double param1, param2, param3, param4,param5,param6;
  int i,j;


  if (Mix_type == 0) pot_parameters(file_echoinput);

  if (Type_attr != NONE){
      for (i=0; i<Ncomp; i++){
         for (j=0; j<Ncomp; j++){
             pairPotparams_switch(Type_pairPot,FLUID_FLUID,i,j,&param1,&param2,&param3,&param4,&param5,&param6);
             Rmin_ff[i][j]=pairPot_find_rmin(i,j,param1,param2,param3,param4,param5,param6,Type_pairPot);
     
         }
      }

      for (i=0; i<Ncomp; i++){
         for (j=0; j<Ncomp; j++){
             pairPotparams_switch(Type_pairPot,FLUID_FLUID,i,j,&param1,&param2,&param3,&param4,&param5,&param6);
             Rzero_ff[i][j]=pairPot_find_r_ZeroCut(i,j,param1,param2,param3,param4,param5,param6,Type_pairPot);
         }
      }

     if (Iwrite_screen==SCREEN_VERBOSE && Proc==0) {
        printf("********************************************************************************\n");
        for (i=0; i<Ncomp; i++){
          for (j=0; j<Ncomp; j++){
             pairPotparams_switch(Type_pairPot,FLUID_FLUID,i,j,&param1,&param2,&param3,&param4,&param5,&param6);
             printf("Rmin_ff[i=%d][j=%d](sigma units)=%15.11f:: check Rmin_ij/sig_ij=%15.11f (const for LJ)\n",i,j,Rmin_ff[i][j],Rmin_ff[i][j]/param1);
          }
       }
        printf("\n");
        for (i=0; i<Ncomp; i++){
          for (j=0; j<Ncomp; j++){
             pairPotparams_switch(Type_pairPot,FLUID_FLUID,i,j,&param1,&param2,&param3,&param4,&param5,&param6);
             printf("Rzero_ff[i=%d][j=%d](sigma units)=%15.11f:: check Rzero_ij/sig_ij=%15.11f (const for LJ if rcut=(const)*sig_ij)\n",i,j,Rzero_ff[i][j],Rzero_ff[i][j]/param1);
          }
       }
       printf("********************************************************************************\n");
     }
  }

  if ((Iwrite_files==FILES_EXTENDED||Iwrite_files==FILES_DEBUG) && Proc==0){
     if (Type_attr != NONE){
       for (i=0; i<Ncomp; i++){
         for (j=0; j<Ncomp; j++){
             print_potentials_fluid(Type_pairPot,i,j);
     
         }
       }
     }
  }
  return;
}
/*******************************************************************************************/
/*pairPot_find_rmin....This function locates the distance r where the minimum in the potential is found.
Note that the distances here are returned in units of r/sigma_ref so no adjustment is needed later.*/

double pairPot_find_rmin(int i, int j,double param1, double param2, double param3,double param4,double param5,double param6,int typePairPot){

  int logical_update=FALSE,count;
  double r,rcut,rmin,umin,umin_old,rstart,delr,sigma,tol=1.e-8,error=1.0,u,uleft,uright,delr_old,rlast;

  rcut=param3;
  sigma=param1;

                                               /* first find general vicinity of minimum */ 
  r=2.0*Esize_x[0];
  delr=Esize_x[0];
  umin=pairPot_switch(r, param1,param2,param3,param4,param5,param6,typePairPot);
  while (r <rcut){
      r+=delr;
      u=pairPot_switch(r, param1,param2,param3,param4,param5,param6,typePairPot);
      if (u<umin){ 
          umin=u;
          rmin=r;
          logical_update=TRUE;
      }
      rlast=r;
  }

  if (logical_update==FALSE || fabs(rmin-rlast)<1.e-6){  
                                 /* no minimum found - this is potential is monotonically increasing or decreasing.  
                                    In either case, the mean field potential term will be defined using a minimum distance of Sigma_ff */
      return(sigma);
  }
  else{                          /* now refine the position of r_min */
      count=1;

      while(error>tol){
           umin_old=umin;
           delr_old=delr; 
           uleft=pairPot_switch(rmin-delr_old, param1,param2,param3,param4,param5,param6,typePairPot);
           uright=pairPot_switch(rmin+delr_old, param1,param2,param3,param4,param5,param6,typePairPot);
       
           rstart=rmin-delr_old;
           r=rstart;
           umin=uleft;

           delr=delr_old/100; 
           while (r<rstart+2.0*delr_old){
               r+=delr;
               u=pairPot_switch(r, param1,param2,param3,param4,param5,param6,typePairPot);
               if (u<umin){ 
                  umin=u;
                  rmin=r;
               }
           }
           error=100*fabs(umin-umin_old)/fabs(umin);
      }
  }
  return(rmin);
}

/*******************************************************************************************/

/*pairPot_find_rzero_cut....This function locates the distances where the interaction potential is equal to the
  interaction potential at r=r_cut.  In these cases u_att=0 for r<r_zero_cut 
  Note that the distances here are returned in units of r/sigma_ref so no adjustment is needed later.*/
double pairPot_find_r_ZeroCut(int i, int j,double param1, double param2, double param3,double param4,double param5,double param6,int typePairPot){

  int logical_update=FALSE,count;
  double r,rcut,rmin,umin,umin_old,rstart,delr,sigma,tol=1.e-8,error=1.0,u,uleft,uright,delr_old,rlast;

  rcut=param3;
  sigma=param1;

                                               /* first find general vicinity of minimum */ 
  r=2.0*Esize_x[0];
  delr=Esize_x[0];
  umin=pairPot_switch(r, param1,param2,param3,param4,param5,param6,typePairPot);
  while (r <rcut-delr){
      r+=delr;
      u=pairPot_switch(r, param1,param2,param3,param4,param5,param6,typePairPot);
      if (u*umin<0.0){ 
          umin=u;
          rmin=r;
          logical_update=TRUE;
      }
      rlast=r;
  }

  if (logical_update==FALSE){  
                                 /* cut and shifted potential never crossed zero so the potential is monotonically increasing or decreasing.  
                                    In either case, the mean field potential term will be defined using a minimum distance of Sigma_ff */
      return(sigma);
  }
  else{                          /* now refine the position of r_min */
      count=1;

      while(error>tol){
           umin_old=umin;
           delr_old=delr; 
           uleft=pairPot_switch(rmin-delr_old, param1,param2,param3,param4,param5,param6,typePairPot);
           uright=pairPot_switch(rmin+delr_old, param1,param2,param3,param4,param5,param6,typePairPot);
       
           rstart=rmin-delr_old; 
           umin=uleft;
           r=rstart;

           delr=delr_old/100; 
           while (r<rstart+2.0*delr_old){
               r+=delr;
               u=pairPot_switch(r, param1,param2,param3,param4,param5,param6,typePairPot);
               if (u*umin<0.0 && fabs(umin)>tol){ 
                  umin=u;
                  rmin=r;
               }
           }
           error=fabs(umin);
      }
  }
  return(rmin);
}
/*******************************************************************************************/
/* print_potentials_fluid  .... output potentials so they can be tested or fitted */
void print_potentials_fluid(int type_pairPot,int icomp,int jcomp){

  double param1, param2, param3, param4,param5,param6,r,uatt,ucs,delr;
  char filename[FILENAME_LENGTH], filenameATT[FILENAME_LENGTH], filenameCS[FILENAME_LENGTH];
  FILE *fpATT, *fpCS;

  pairPotparams_switch(Type_pairPot,FLUID_FLUID,icomp,jcomp,&param1,&param2,&param3,&param4,&param5,&param6);

  sprintf(filenameATT, "dft_uATT%0d%0d.dat", icomp,jcomp);
  sprintf(filenameCS, "dft_uCS%0d%0d.dat", icomp,jcomp);
  fpATT=fopen(filenameATT,"w");
  fpCS=fopen(filenameCS,"w");
 
  delr=Esize_x[0]/10; 
  for (r=delr; r<Cut_ff[icomp][jcomp]+1.0; r+=delr){
     uatt=pairPot_ATT_CS_switch(r,icomp,jcomp,type_pairPot);
     fprintf(fpATT,"%11.6f  %11.6f\n",r,uatt); 

     ucs=pairPot_switch(r,param1,param2,param3,param4,param5,param6,type_pairPot);
     fprintf(fpCS,"%11.6f  %11.6f\n",r,ucs); 
  }
  return; 
}
/*******************************************************************************************/
