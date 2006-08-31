/*
//@HEADER
// ********************************************************************
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER
*/

/* ---------------------------------------------------------
Calculate relevant thermodynamic properties for CMS polymer fluid.
------------------------------------------------------------*/
#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
/***************************************************************
/*setup_polymer_cr: read in c(r) from file and add attractions*/
void setup_polymer_cr()
{     
   FILE *fp7,*fp8,*fp9,*fp10;
   int i,j,lines,lines2,ir,ir2, nextChar;
   double r,u,cr_rad_max,cr_rad_max2,rsave,rsave2,dummy_read,crread;
   double crfac1,crfac2,crfac3,crfac4,xs=0.0;
   char c;

   /* note that this routine is set up to use either a single direct correlation function as input - or
      to automatically interpolate between up to 4 files.  The particular automated interpolation set up
      now involves the fraction of solvent in the bulk uniform fluids, xs.  This was useful in the
      studies of lipid bilayers (Frink and Frischknecht, PRE 2005) Other interpolations would require
      editing this file */
     
      if (Ncomp==3){ 
        xs=Rho_b[2]/(Rho_b[0]+Rho_b[1]+Rho_b[2]);
      }
      else{
        if (Ncr_files>1){ 
             printf("Number of cr_files must be 1 unless doing a 3 component system at this point\n");
             exit(-1);
         }
      } 
      if (Ncr_files == 1){
             crfac1=Crfac;
      }
      else if (Ncr_files == 2){
             crfac1=xs;
             crfac2=1.0-xs;
      }
      else if (Ncr_files ==3){
         if (xs>0 && xs<Cr_break[0]){
            crfac1=xs/(Cr_break[0]);
            crfac2=1.0-crfac1;
            crfac3=0.0;
         }
         else{
           crfac1=0.0;
           crfac2=(xs-Cr_break[0])/(1.0-Cr_break[0]);
           crfac3=1.0-crfac2;
         }
      }
      else if (Ncr_files ==4){
         if (xs>0 && xs<Cr_break[0]){
            crfac1=1.0-xs/(Cr_break[0]);
            crfac2=1.0-crfac1;
            crfac3=0.0; crfac4=0.0;
         }
         else if (xs>=Cr_break[0] && xs<Cr_break[1]){
            crfac1=0.0; crfac4=0.0;
            crfac2 = 1.0-(xs-Cr_break[0])/(Cr_break[1]-Cr_break[0]);
            crfac3=1.0-crfac2;
         }
         else if (xs>=Cr_break[1]){
            crfac1=0.0; crfac2=0.0;
            crfac3 = 1.0 - (xs-Cr_break[1])/(1.0-Cr_break[1]);
            crfac4=1.0-crfac3;
         }
      }
   if (Proc==0 && Iwrite==VERBOSE){
      if (Ncr_files >= 1) printf("crfac1=%9.6f  ",crfac1);
      if (Ncr_files >= 2) printf("crfac2=%9.6f  ",crfac2);
      if (Ncr_files >= 3) printf("crfac3=%9.6f  ",crfac3);
      if (Ncr_files >= 4) printf("crfac4=%9.6f  ",crfac4);
      printf("\n");
   }

   /* reading in c(r) file */
   if(Proc==0) printf("reading in %d c(r) file(s)...\n",Ncr_files);
   if (Type_poly == CMS) {

   if (Proc==0){
     if( (fp7  = fopen(Cr_file,"r")) == NULL) {
       printf("Can't open file %s\n", Cr_file);
       exit(1);
     }
     if (Ncr_files>=2){
       if( (fp8  = fopen(Cr_file2,"r")) == NULL) {
          printf("Can't open file %s\n", Cr_file2);
          exit(1);
        }
        fclose(fp8);
     }
     if (Ncr_files>=3){
       if( (fp9  = fopen(Cr_file3,"r")) == NULL) {
          printf("Can't open file %s\n", Cr_file3);
          exit(1);
        }
        fclose(fp9);
     }
     if (Ncr_files==4){
       if( (fp10  = fopen(Cr_file4,"r")) == NULL) {
          printf("Can't open file %s\n", Cr_file4);
          exit(1);
        }
        fclose(fp10);
     }
     for (ir=1; ir<=3; ir++){
       fscanf(fp7,"%lf",&r );
       fscanf(fp7,"%c",&c );
       for (i=0; i<Ncomp; i++)  /* for (i=0; i<Ntype_mer; i++)  */
         fscanf(fp7,"%lf",&dummy_read);
       fscanf(fp7,"%c",&c );
       for (i=0; i<Ncomp; i++) {                /* "  */
         for (j=i+1; j<Ncomp; j++) {                /* "  */
           fscanf(fp7,"%lf",&dummy_read);
           fscanf(fp7,"%c",&c );
         }
       }
       if (ir == 1) rsave = r;
       if (ir == 2) Deltar_cr = r-rsave;
     }
     fclose(fp7);
   }
   MPI_Bcast(&Deltar_cr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   lines = 0;
   cr_rad_max=0.0;
   for (i=0; i < Ncomp; ++i)
      for (j=0; j<Ncomp; ++j){
         if (Cr_rad_hs[i][j] > cr_rad_max) cr_rad_max = Cr_rad_hs[i][j];
      }
   lines=(int)(cr_rad_max/Deltar_cr);

   if (Proc==0 && lines >= N_NZCR_MAX) {
      printf("Warning: Rism_cr array may be truncated :: Max Cr_rad_hs > allowed\n");
      lines = N_NZCR_MAX-1;
    }
    /*printf("Proc=%d lines=%d\n",Proc,lines); */

   Last_nz_cr = lines;
   if (Proc==0) {
     fp7  = fopen(Cr_file,"r");
     if (Ncr_files>=2) fp8  = fopen(Cr_file2,"r");
     if (Ncr_files>=3) fp9  = fopen(Cr_file3,"r");
     if (Ncr_files==4) fp10  = fopen(Cr_file4,"r");

     for (ir=1; ir<=lines; ir++){
       fscanf(fp7,"%lf",&r );
       fscanf(fp7,"%c",&c );
       if (ir == 1) rsave = r;
       if (ir == 2) Deltar_cr = r-rsave;

       if (Ncr_files>=2){
         fscanf(fp8,"%lf",&r );
         fscanf(fp8,"%c",&c );
       }
       if (Ncr_files>=3){
         fscanf(fp9,"%lf",&r );
         fscanf(fp9,"%c",&c );
       }
       if (Ncr_files==4){
         fscanf(fp10,"%lf",&r );
         fscanf(fp10,"%c",&c );
       }

       for (i=0; i<Ncomp; i++){  /* for (i=0; i<Ntype_mer; i++)  */
         fscanf(fp7,"%lf",&crread);
         Rism_cr[i][i][ir]=Crfac*crfac1*crread;
         if (Ncr_files>=2){
              fscanf(fp8,"%lf",&crread);
              Rism_cr[i][i][ir]+=Crfac*crfac2*crread;
         }
         if (Ncr_files>=3){
              fscanf(fp9,"%lf",&crread);
              Rism_cr[i][i][ir]+=Crfac*crfac3*crread;
         }
         if (Ncr_files==4){
              fscanf(fp10,"%lf",&crread);
              Rism_cr[i][i][ir]+=Crfac*crfac4*crread;
         }
       }
       fscanf(fp7,"%c",&c );
       if (Ncr_files>=2) fscanf(fp8,"%c",&c );
       if (Ncr_files>=3) fscanf(fp9,"%c",&c );
       if (Ncr_files==4) fscanf(fp10,"%c",&c );

       for (i=0; i<Ncomp; i++) {
         for (j=i+1; j<Ncomp; j++) {
           fscanf(fp7,"%lf",&crread);
           Rism_cr[i][j][ir]=Crfac*crfac1*crread;
           fscanf(fp7,"%c",&c );
           if (Ncr_files>=2){
                 fscanf(fp8,"%lf",&crread);
                 Rism_cr[i][j][ir]+=Crfac*crfac2*crread;
                 fscanf(fp8,"%c",&c );
           }
           if (Ncr_files>=3){
                 fscanf(fp9,"%lf",&crread);
                 Rism_cr[i][j][ir]+=Crfac*crfac3*crread;
                 fscanf(fp9,"%c",&c );
           }
           if (Ncr_files==4){
                 fscanf(fp10,"%lf",&crread);
                 Rism_cr[i][j][ir]+=Crfac*crfac4*crread;
                 fscanf(fp10,"%c",&c );
           }
           Rism_cr[j][i][ir] = Rism_cr[i][j][ir];
         }
       }
       while(c != '\n') c=getc(fp7);
       if(Ncr_files>=2) while(c != '\n') c=getc(fp8);
       if(Ncr_files>=3) while(c != '\n') c=getc(fp9);
       if(Ncr_files==4) while(c != '\n') c=getc(fp10);
     }
     fclose(fp7);
     if(Ncr_files>=2) fclose(fp8);
     if(Ncr_files>=3) fclose(fp9);
     if(Ncr_files==4) fclose(fp10);

     for (ir=lines+1; ir<N_NZCR_MAX; ir++)
       for (i=0; i<Ncomp; i++){
         Rism_cr[i][i][ir] = 0.;
         for (j=i+1; j<Ncomp; j++) {
           Rism_cr[i][j][ir] = 0.;
           Rism_cr[j][i][ir] = 0.;
         }
       }
     /* extrapolate c(r) to r = 0 */
     for (i=0; i<Ncomp; i++)                 /* "  */
       for (j=0; j<Ncomp; j++)                 /* "  */
         Rism_cr[i][j][0] = 2.*Rism_cr[i][j][1] - Rism_cr[i][j][2];
   }
   MPI_Bcast(**Rism_cr,Ncomp*Ncomp*N_NZCR_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /* add attractions to polymer c(r)  */

   if (Iwrite == VERBOSE && Proc==0){
      fp7  = fopen("cr.out","w");
      for (ir=0; ir <= Last_nz_cr; ++ir){
         fprintf(fp7,"%lf   ",ir*Deltar_cr);
         for (i=0; i < Ncomp; ++i)
            for (j=0; j < Ncomp; ++j)
               fprintf(fp7,"%lf   ",Rism_cr[i][j][ir]);
         fprintf(fp7,"\n");
      }
      fclose(fp7);
   }

   if (Ipot_ff_n == LJ12_6){
     for (i=0; i < Ncomp; ++i)
       for (j=0; j<Ncomp; ++j){
         if (Cut_ff[i][j] > Cr_rad_hs[i][j]) Cr_rad[i][j] = Cut_ff[i][j];
         else                                Cr_rad[i][j] = Cr_rad_hs[i][j];
         if (Cr_rad[i][j] > cr_rad_max) cr_rad_max=Cr_rad[i][j];
       }

   lines = cr_rad_max/Deltar_cr;
   if (lines >= N_NZCR_MAX) {
     if (Proc==0) printf("Need to increase N_NZCR_MAX\n");
     lines = N_NZCR_MAX-1;
   }
   if (lines > Last_nz_cr) Last_nz_cr = lines;

     for (i=0; i < Ncomp; ++i)
       for (j=i; j<Ncomp; j++) {
         for (ir=0; ir<=lines; ir++){
           r = ir*Deltar_cr;
           if (r-1.e-8 > Sigma_ff[i][j]*1.122462) { /* watch roundoffs*/
             u = uLJatt_n(r,i,j);
             Rism_cr[i][j][ir] -= u;
             if (i != j) Rism_cr[j][i][ir] -= u;
           }
         }
       }
   }
   if (Iwrite == VERBOSE && Proc==0){
      fp7  = fopen("cr.lj.out","w");
      for (ir=0; ir <= Last_nz_cr; ++ir){
         fprintf(fp7,"%lf   ",ir*Deltar_cr);
         for (i=0; i < Ncomp; ++i)
            for (j=0; j < Ncomp; ++j)
               fprintf(fp7,"%lf   ",Rism_cr[i][j][ir]);
         fprintf(fp7,"\n");
      }
      fclose(fp7);
   }

   }
   else if (Type_poly==CMS_SCFT) { /* here only read c(0)*/
     if (Proc==0){
       if( (fp7  = fopen(Cr_file,"r")) == NULL) {
         printf("Can't open file %s\n", Cr_file);
         exit(1);
       }

       for (i=0; i<Ncomp; i++) {  /* for (i=0; i<Ntype_mer; i++)  */
         fscanf(fp7,"%lf",&Rism_cr[i][i][0]);
         fscanf(fp7,"%c",&c );
         printf("cr[%d][%d] = %f\t", i, i, Rism_cr[i][i][0]);
       }
       for (i=0; i<Ncomp; i++) {                /* "  */
         for (j=i+1; j<Ncomp; j++) {                /* "  */
           fscanf(fp7,"%lf",&Rism_cr[i][j][0]);
           fscanf(fp7,"%c",&c );
           Rism_cr[j][i][0] = Rism_cr[i][j][0];
           printf("cr[%d][%d] = %f\t", i,j,Rism_cr[i][j][0]);
         }
       }
       fclose(fp7);
       for (ir=1; ir<N_NZCR_MAX; ir++){
         for (i=0; i<Ncomp; i++){
           Rism_cr[i][i][ir] = 0.;
           for (j=i+1; j<Ncomp; j++) {
             Rism_cr[i][j][ir] = 0.;
             Rism_cr[j][i][ir] = 0.;
           }
         }
       }
     } /* end of if(Proc==0) */

     MPI_Bcast(**Rism_cr,Ncomp*Ncomp*N_NZCR_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   return;
}

