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
 *  FILE: dft_input.c
 *
 *  This file contains routines that read in the input file and
 *  initializes some flags based on the input.
 */


#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "dft_globals_const.h"
#include "rf_allo.h"
#include <mpi.h>

/* Prototypes for functions found in this file */
void read_junk(FILE *fp, FILE *fp2);  
void error_check(void);

/************** R O U T I N E S   I N   T H E   F I L E ***********************
*
*    NAME                               TYPE            CALLED_BY
*--------------------------------------------------------------------
*
* read_input_file                       void            main (dft_main.c)
* read_junk				void		read_input_file
*
******************************************************************************/

void read_input_file(char *input_file, char *output_file1)

/* Read the input file for Classical Fluids Density Functional Theory
 *
 *    Authors:  Laura Frink, 9225
 *              Andrew Salinger, 9221
 */

{
   /* Local variable declarations */
   
   FILE *fp=NULL, *fp2=NULL, *fp3=NULL, *fp4=NULL;

   char *yo = "read_input_file";
   char poly_file[20];
   int isten, icomp, jcomp, iwall,iwall_save,iwall_type, idim, 
       i, izone, j, jwall,end_count,
       new_wall,logical,ncharge, seg, block[NCOMP_MAX][NBLOCK_MAX],
       block_type[NBLOCK_MAX],pol_number, nlink_chk,irand,irand_range,itmp,
       *nbond_tot,nbond_all,iseg,nseg,nmer_max,ibond,pol_num2,nunk,
       ***pol_sym_tmp,dim_tmp,Lauto_center,Lauto_size,Lcompare_fastram,jmin=0,jmax=0,
       lzeros,latoms,ltrues,jwall_type,seg_tot;
   double r,rho_tmp[NCOMP_MAX],dxdx,dtmp,charge_sum,minpos[3],maxpos[3];
   double rough_param_max[NWALL_MAX_TYPE],rough_length_scale[NWALL_MAX_TYPE];
   int iblock,jblock;

  
  /********************** BEGIN EXECUTION ************************************/

  if (Proc==0) {
     printf("\n-----------------STARTING DFT CALCULATION----------------------------\n");
  }
  /* Read in the Mesh, Surface, Potential Type, Fluid Particle, 
     Surface Particle, State Point, Functional, and Run Control Parameters.

     The Input file (dft_input.dat) may be freely formatted with comments
     as long as the @ symbol is placed at the beginning of all data input
     lines.

     The input file is copied to the output file dft_out.lis */

  /* Open the Files ALF also check for errors! */
  if (Proc==0) {
    if( (fp  = fopen(input_file,"r")) == NULL) {
      printf("Can't open file %s\n", input_file);
      exit(1);
    }
    if( (fp2 = fopen(output_file1,"w+")) == NULL) {
      printf("Can't open file %s\n", output_file1);
      exit(1);
    }
    fprintf(fp2,"test out the printing %s\n",output_file1);
  }

  /* Initialize and Read Dimension Parameters */
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%lf  %lf  %lf  %lf %lf", 
             &Length_ref, &Density_ref, &Temp, &Dielec_ref, &VEXT_MAX);
    fprintf(fp2,"%f  %f  %f  %f   %f",
               Length_ref,Density_ref,Temp,Dielec_ref,VEXT_MAX);
  }
  MPI_Bcast(&Length_ref,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Density_ref,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Temp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Dielec_ref,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&VEXT_MAX,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Temp < 0.0) Potential_ref=-1.0;
  else            Potential_ref=1000.*KBOLTZ*Temp/E_CONST;  /* reference potential in mV */

  /* Initialize and Read Mesh Parameters */

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Ndim);
    fprintf(fp2,"%d",Ndim);
  }
  MPI_Bcast(&Ndim,1,MPI_INT,0,MPI_COMM_WORLD);


  if (Proc==0) {
    read_junk(fp,fp2);
    for (idim=0; idim < Ndim; ++idim){
      fscanf(fp,"%lf", &Size_x[idim]);
      fprintf(fp2,"%f  ",Size_x[idim]);
      if (Length_ref >0.0) Size_x[idim] /= Length_ref;
    }
    for (idim=Ndim; idim<NDIM_MAX; idim++) Size_x[idim]=0.0;
  }
  MPI_Bcast(Size_x,NDIM_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    for (idim=0; idim < Ndim; ++idim){
      fscanf(fp,"%lf", &Esize_x[idim]);
      fprintf(fp2,"%f  ",Esize_x[idim]);
      if (Length_ref > 0.0) Esize_x[idim] /= Length_ref;
    }
    for (idim=Ndim; idim<NDIM_MAX; idim++) Esize_x[idim]=0.0;
  }
  MPI_Bcast(Esize_x,NDIM_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

/*  read_junk(fp,fp2);
 * fscanf(fp,"%d",&Lmesh_refine);
 * if (Proc==0) fprintf(fp2,"%d",Lmesh_refine);
 */ Lmesh_refine=FALSE; /*parameter not active currently - auto mesh refine */

  if ( Proc==0 ) {
    for (idim=0; idim < Ndim; ++idim){
      read_junk(fp,fp2);
      for (i=0; i < 2; i++){
	fscanf(fp,"%d", &Type_bc[idim][i]);
	fprintf(fp2,"%d  ",Type_bc[idim][i]);
      }
    }
  }
  MPI_Bcast(Type_bc,Ndim*2,MPI_INT,0,MPI_COMM_WORLD);
  if ( Proc==0 ) {
    if (Ndim == 1) {
      read_junk(fp,fp2);
      read_junk(fp,fp2);
    }
    else if (Ndim ==2) {
      read_junk(fp,fp2);
    }
  }
  for (idim=0; idim < Ndim; ++idim) {
    if ((Type_bc[idim][0]==PERIODIC)-(Type_bc[idim][1]==PERIODIC)) {
         printf("%s: ERROR: Both BC in dimension %d must be periodic if one is\n",
              yo, idim);
         exit(-1); 
    }
  }


  /* Read in Functional Switches  and set up the Sten_Type array.*/
  for (isten=0; isten<NSTEN; ++isten) Sten_Type[isten]=FALSE;

  /* hard sphere functionals */
  if ( Proc==0 ) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Type_func);
    fprintf(fp2,"%d",Type_func);
  }
  MPI_Bcast(&Type_func,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Type_func >=0 && Type_func <= 2){
    Sten_Type[DELTA_FN]=Sten_Type[THETA_FN]=TRUE;
  }
  else if (Type_func >2 || Type_func<-1){
    if (Proc==0) printf("ERROR Type_hs out of range - should be -1,0,1\n");
    exit(-1);
  }

  /* attractive functionals */
  if ( Proc==0 ) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Type_attr);
    fprintf(fp2,"%d",Type_attr);
  }
  MPI_Bcast(&Type_attr,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Type_attr==0 || Type_attr==1) Sten_Type[U_ATTRACT]=TRUE;
  else if (Type_attr >1 || Type_attr<-1){
     if (Proc==0) printf("ERROR Type_attr=%d out of range - should be -1, 0, or 1\n",Type_attr);
     exit(-1);
  }

  /* coulombic functionals */
  if ( Proc==0 ) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Type_coul);
    fprintf(fp2,"%d",Type_coul);
  }
  MPI_Bcast(&Type_coul,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Type_coul==1) Sten_Type[THETA_CHARGE]=TRUE;
  else if (Type_coul >2 || Type_coul<-1){
    if (Proc==0) printf("ERROR Type_coul out of range - should be -1,0,1\n");
    exit(-1);
  }

  /* polymer functionals */
  if ( Proc==0 ) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Type_poly);
    fprintf(fp2,"%d",Type_poly);
  }
  MPI_Bcast(&Type_poly,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Type_poly == WTC) {
    Sten_Type[THETA_FN_SIG]=TRUE;
    Sten_Type[DELTA_FN_BOND]=TRUE;
  }
  else if (Type_poly == CMS || Type_poly==CMS_SCFT){
      Sten_Type[DELTA_FN]=TRUE;
      Sten_Type[U_ATTRACT]=FALSE;   /* attractions handled differently for polymers */
      Sten_Type[POLYMER_CR]=TRUE;
      if (Type_poly==CMS_SCFT){
        printf ("To do SCFT with CMS theory, we need to test and debug all code !\n");
        exit(-1);
      }
  }
  if (Type_poly == CMS_GAUSSIAN){
      Sten_Type[POLYMER_CR]=2;
      printf ("To do CMS Gaussian chains, we need to test and debug all code !\n");
      exit(-1);
  }
  else if (Type_poly >WTC || Type_poly<NONE){
     if (Proc==0) printf("ERROR Type_poly out of range (bounds are %d,%d)\n",NONE,WTC);
     exit(-1);
  }

/* check if we are trying to compare with FasTram -- for hard systems */
  if ( Proc==0 ) {
    read_junk(fp,fp2);
    fscanf(fp,"%d ",&Lcompare_fastram);
    fprintf(fp2,"%d ",Lcompare_fastram);
  }
  MPI_Bcast(&Lcompare_fastram,1,MPI_INT,0,MPI_COMM_WORLD);

  /* Read in or set if known the Potential Type Paramters */
  if (Type_func == -1 && (Type_poly==NONE || Type_poly==WTC)  ) Ipot_ff_n = IDEAL_GAS;
  else if (Type_attr == -1)                                     Ipot_ff_n = HARD_SPHERE;
  else if (Type_attr == 0 || Type_attr==1)                                      Ipot_ff_n = LJ12_6;
  else {
     printf("ERROR WITH Type_func and Type_attr selections and conversion to Ipot_ff_n parameter \n");
     exit (-1);
  }
  MPI_Bcast(&Ipot_ff_n,1,MPI_INT,0,MPI_COMM_WORLD);

 
  if (Type_coul >=0)  Ipot_ff_c =1; /* Coulombic Fluid */
  else Ipot_ff_c=0;  /* Neutral Fluid */
  MPI_Bcast(&Ipot_ff_c,1,MPI_INT,0,MPI_COMM_WORLD);

  /* Read in Surface Paramters */
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d  %d  %d  %d  %d",&Nwall_type,&Nwall,&Nlink,&Lauto_center,&Lauto_size);
    fprintf(fp2,"%d  %d  %d  %d  %d",Nwall_type,Nwall,Nlink,Lauto_center,Lauto_size);
  }
  MPI_Bcast(&Nwall_type,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Nwall,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Nlink,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Lauto_center,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Lauto_size,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Nwall>0) Xtest_reflect_TF = (int **) array_alloc (2, Nlink,Ndim, sizeof(int));
  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall > 0)
      for (i=0; i < Nlink; i++)
	for (idim=0; idim< Ndim; idim++){
	  fscanf(fp,"%d",&Xtest_reflect_TF[i][idim]);
	  fprintf(fp2,"%d  ",Xtest_reflect_TF[i][idim]);
	}
  }
  if (Nwall>0)
  MPI_Bcast(*Xtest_reflect_TF,Nlink*Ndim,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%d", &Surface_type[iwall_type]);
	fprintf(fp2,"%d  ",Surface_type[iwall_type]);
      }
    else fprintf(fp2,"n/a");
  }
  if(Nwall_type > 0) 
    MPI_Bcast(Surface_type,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%d", &Orientation[iwall_type]);
	fprintf(fp2,"%d  ",Orientation[iwall_type]);
      }
    else fprintf(fp2,"n/a");
  }
  if (Nwall_type > 0) 
    MPI_Bcast(Orientation,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);


  charge_sum=0.0;
  srandom(135649);
  if (Nwall_type != 0){
    if ( Proc==0) {
      for (idim=0; idim<Ndim; idim++){ minpos[idim] = 1000.; maxpos[idim]=-1000.;}
      /* ALF: add error checking */
      if( (fp3  = fopen("dft_surfaces.dat","r")) == NULL) {
	printf("Can't open file dft_surfaces.dat\n");
	exit(1);
      }
      for (iwall=0; iwall<Nwall; iwall++){
	fscanf(fp3,"%d  %d",&WallType[iwall], &Link[iwall]);
	for (idim=0; idim<Ndim; idim++) {
          dim_tmp=idim;
  /* temporary rotation of coordinates 
     if (idim==0) dim_tmp=1;
     else if (idim==1) dim_tmp=2;
     else if (idim==2) dim_tmp=0;*/
  /* end of temporary code */
	  fscanf(fp3,"%lf",&WallPos[dim_tmp][iwall]);
          if (fabs(WallPos[dim_tmp][iwall]+9999.0)<1.e-6) { /*random coordinate placement */
             irand = random();
             irand_range = POW_INT(2,31)-1;
              WallPos[dim_tmp][iwall] = Size_x[idim]*(-0.5+( ((double)irand)/((double)irand_range)));
              printf("\n  Wall %d dim %d gets WallPos:%g \n",iwall,idim,WallPos[idim][iwall]);
              fprintf(fp2,"\n Wall %d dim %d gets WallPos:%g \n",iwall,idim,WallPos[idim][iwall]);
          }
          if (Length_ref > 0.0) WallPos[dim_tmp][iwall] /=Length_ref; 
          if (WallPos[dim_tmp][iwall] < minpos[dim_tmp]) minpos[dim_tmp]=WallPos[dim_tmp][iwall];
          if (WallPos[dim_tmp][iwall] > maxpos[dim_tmp]) maxpos[dim_tmp]=WallPos[dim_tmp][iwall];
       }
       fscanf(fp3,"%lf",&Elec_param_w[iwall]);
       charge_sum+=Elec_param_w[iwall];
       } 
       /*for (idim=0; idim<Ndim; idim++) printf("\n idim: %d min pos: %9.6f max pos %9.6f \n",idim,minpos[idim],maxpos[idim]);*/
       if (Lauto_center){
         for (iwall=0; iwall<Nwall; iwall++)
            for (idim=0; idim<Ndim; idim++) WallPos[idim][iwall] -= 0.5*(maxpos[idim] + minpos[idim]);
       }

      fclose(fp3);
      }
    MPI_Bcast(WallType,NWALL_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Link,NWALL_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(WallPos,NDIM_MAX*NWALL_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if (Nwall_type>0)
    MPI_Bcast(Elec_param_w,NWALL_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

    nlink_chk = 1;
    for (iwall=1; iwall<Nwall; iwall++){
      new_wall = TRUE;
      for (jwall=0; jwall<iwall; jwall++)
	if (Link[iwall] == Link[jwall]) new_wall=FALSE;
      if (new_wall) nlink_chk++;
    }
    if (nlink_chk != Nlink){
      printf("Check Nlink in dft_input.dat: %d and assignments in dft_surfaces.dat: %d\n",
	     Nlink,nlink_chk);
      exit(-1);
    }
    Link_list = (int **) array_alloc (2, Nlink,Nwall,sizeof(int));
    Nwall_this_link = (int *) array_alloc (1, Nlink,sizeof(int));
    for (i=0; i<Nlink; i++) Nwall_this_link[i]=0;
    for (iwall=0; iwall<Nwall; iwall++){
      Link_list[Link[iwall]][Nwall_this_link[Link[iwall]]++]=iwall;
    }
  }
  else{
    Nwall = 0;
    Nlink = 0;
  }

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%lf", &WallParam[iwall_type]);
	fprintf(fp2,"%f  ",WallParam[iwall_type]);
        if (Length_ref >0.0) WallParam[iwall_type]/=Length_ref;
        if (Surface_type[iwall_type]==point_surface) WallParam[iwall_type]=Esize_x[0];
      }
    else fprintf(fp2,"n/a");
  }
  if (Nwall_type > 0) 
    MPI_Bcast(WallParam,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%lf", &WallParam_2[iwall_type]);
	fprintf(fp2,"%f  ",WallParam_2[iwall_type]);
        if (Length_ref >0.0) WallParam_2[iwall_type]/=Length_ref;
      }
    else fprintf(fp2,"n/a");
  }
  
  if (Nwall_type > 0) 
     MPI_Bcast(WallParam_2,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%lf", &WallParam_3[iwall_type]);
	fprintf(fp2,"%f  ",WallParam_3[iwall_type]);
        if (Length_ref >0.0) WallParam_3[iwall_type]/=Length_ref;
      }
    else fprintf(fp2,"n/a");
  }
  if (Nwall_type > 0) 
    MPI_Bcast(WallParam_3,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%lf", &WallParam_4[iwall_type]);
	fprintf(fp2,"%f  ",WallParam_4[iwall_type]);
        if (Length_ref >0.0) WallParam_4[iwall_type]/=Length_ref;
      }
    else fprintf(fp2,"n/a");
  }
  if (Nwall_type > 0) 
    MPI_Bcast(WallParam_4,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);


  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%d", &Lrough_surf[iwall_type]);
	fprintf(fp2,"%d  ",Lrough_surf[iwall_type]);
      }
    else fprintf(fp2,"n/a");
  }
  if (Nwall_type > 0) 
    MPI_Bcast(Lrough_surf,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%lf", &rough_param_max[iwall_type]);
	fprintf(fp2,"%f  ",rough_param_max[iwall_type]);
      }
    else fprintf(fp2,"n/a");
  }
  if (Nwall_type > 0) 
    MPI_Bcast(rough_param_max,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%lf", &Rough_length[iwall_type]);
	fprintf(fp2,"%l  ",Rough_length[iwall_type]);
      }
    else fprintf(fp2,"n/a");
  }
  if (Nwall_type > 0) 
    MPI_Bcast(Rough_length,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

         /* we don't really know the correct number of rouch tiles (blocks) 
           for any problem at this point .... so we will just populate the random roughness array fully */
  for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){
      for (iblock=0;iblock<MAX_ROUGH_BLOCK;iblock++){
         for (jblock=0;jblock<MAX_ROUGH_BLOCK;jblock++){
            irand = random();
            irand_range = POW_INT(2,31)-1;
            Rough_precalc[iwall_type][iblock][jblock]= rough_param_max[iwall_type]*(-0.5+( ((double)irand)/((double)irand_range)));
         }
      }
  }

  /* switches for types of wall-fluid and wall-wall interaction parameters */
  Lhard_surf=FALSE;
  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         fscanf(fp,"%d",&Ipot_wf_n[iwall_type]);
         fprintf(fp2,"%d",Ipot_wf_n[iwall_type]);
         if (Ipot_wf_n[iwall_type]==VEXT_HARD) Lhard_surf=TRUE;
          /* set logical indicating if any of the surfaces have hard cores - if so, we
              will need be careful with rosenfeld integrals */
      }
    else fprintf(fp2,"n/a");
  }
  if (Lcompare_fastram) Lhard_surf=FALSE;
  MPI_Bcast(Ipot_wf_n,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Lhard_surf,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    lzeros=FALSE; latoms=FALSE;
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
        for (jwall_type=0; jwall_type < Nwall_type;++jwall_type){
           if (!lzeros && !latoms){
             fscanf(fp,"%d",&Ipot_ww_n[iwall_type][jwall_type]);
             fprintf(fp2,"%d",Ipot_ww_n[iwall_type][jwall_type]);
           }
           if (iwall_type==0 && jwall_type==0){
                if (Ipot_ww_n[iwall_type][jwall_type]==-2) lzeros=TRUE;
                else if (Ipot_ww_n[iwall_type][jwall_type]==-1) latoms=TRUE;
           }
           if (lzeros) Ipot_ww_n[iwall_type][jwall_type]=0;
           else if (latoms) Ipot_ww_n[iwall_type][jwall_type]=1;
        }
      }
    else fprintf(fp2,"n/a");
  }
  MPI_Bcast(Ipot_wf_n,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);



  /* Fluid Particle Parameters */

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d  %d",&Ncomp,&Mix_type);
    fprintf(fp2,"%d  %d",Ncomp,Mix_type);
  }
  MPI_Bcast(&Ncomp,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Mix_type,1,MPI_INT,0,MPI_COMM_WORLD);

/* New code for interaction potential parameters */

                      /* MASS Entries */
  if (Proc==0) {
    read_junk(fp,fp2);
    for (i=0; i<Ncomp; i++){
      fscanf(fp,"%lf",&Mass[i]); 
      fprintf(fp2,"%f  ",Mass[i]);
    }
  }
    MPI_Bcast(Mass,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

                      /* VALENCE Entries */
  if (Proc==0) {
    read_junk(fp,fp2);
    if (Ipot_ff_c >0){
      for (icomp=0; icomp < Ncomp; ++icomp){
        fscanf(fp,"%lf", &Charge_f[icomp]);
        fprintf(fp2,"%f  ",Charge_f[icomp]);
      }
    }
    else{ 
      fprintf(fp2,"n/a");
      for (icomp=0; icomp < Ncomp; ++icomp)Charge_f[icomp]=0.0;
    }
  }
  MPI_Bcast(Charge_f,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

                      /* POLARIZATION Entries */
  if (Proc==0) {
    read_junk(fp,fp2);
    if (Type_coul ==2){
       for (icomp=0; icomp < Ncomp; ++icomp){
          fscanf(fp,"%lf", &Pol[icomp]);
          if (Proc==0) fprintf(fp2,"%f  ",Pol[icomp]);
          if (Pol[icomp]!=0.0) Lpolarize[icomp]=TRUE;
          else                 Lpolarize[icomp]=FALSE;
       }
    }
    else if (Proc==0) fprintf(fp2,"n/a: no polarization");
  }
  MPI_Bcast(Pol,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(Lpolarize,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);

                       /* FLUID-FLUID PARAMS */

                       /* Sigma_ff */
  if (Proc==0) {
    read_junk(fp,fp2);
    if (Ipot_ff_n != IDEAL_GAS){
      for (i=0; i<Ncomp; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
        for (j=jmin; j<jmax; j++){
	  fscanf(fp,"%lf",&Sigma_ff[i][j]);
	  fprintf(fp2,"%f  ",Sigma_ff[i][j]);
          if (Length_ref > 0.0) Sigma_ff[i][j]/=Length_ref;
        }
      }
    }
    else  {
      for (i=0; i<Ncomp; i++){
         if (Mix_type==0) {jmin=i; jmax=i+1;}
         else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
         for(j=jmin;j<jmax;j++) Sigma_ff[i][j] = 0.0;
      }
      fprintf(fp2,"n/a");
    }
  }
  MPI_Bcast(Sigma_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

                       /* Eps_ff */
  if (Proc==0) {
    read_junk(fp,fp2);
    if (Ipot_ff_n != IDEAL_GAS && Ipot_ff_n != HARD_SPHERE) {
      for (i=0; i<Ncomp; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
        for (j=jmin; j<jmax; j++){
	  fscanf(fp,"%lf",&Eps_ff[i][j]);
	  fprintf(fp2,"%f  ",Eps_ff[i][j]);
          if (Temp > 0.0) Eps_ff[i][j]/=Temp;
        }
      }
    }
    else  {
      for (i=0; i<Ncomp; i++){ 
         if (Mix_type==0) {jmin=i; jmax=i+1;}
         else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
         for (j=jmin; j<jmax; j++) Eps_ff[i][j] = 0.0;
      }
      fprintf(fp2,"n/a: no Eps_ff for ideal gas or hard sphere");
    }
  }
  MPI_Bcast(Eps_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

                       /* Cut_ff */
  if (Proc==0) {
    read_junk(fp,fp2);
    if (Ipot_ff_n != IDEAL_GAS && Ipot_ff_n != HARD_SPHERE) {
      for (i=0; i<Ncomp; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
        for (j=jmin; j<jmax; j++){
	  fscanf(fp,"%lf",&Cut_ff[i][j]);
	  fprintf(fp2,"%f  ",Cut_ff[i][j]);
          if (Length_ref > 0.0) Cut_ff[i][j]/=Length_ref;
        }
      }
    }
    else  {
      for (i=0; i<Ncomp; i++){ 
         if (Mix_type==0) {jmin=i; jmax=i+1;}
         else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
         for (j=jmin; j<jmax; j++) Cut_ff[i][j] = 0.0;
      }
      fprintf(fp2,"n/a: no cutoff for ideal gas or hard sphere");
    }
  }
  MPI_Bcast(Cut_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

                       /* Bond_ff */
  if (Proc==0) {
    read_junk(fp,fp2);
    if (Type_poly!=NONE){
      for (i=0; i<Ncomp; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
        for (j=jmin; j<jmax; j++){
	  fscanf(fp,"%lf",&Bond_ff[i][j]);
	  fprintf(fp2,"%f  ",Bond_ff[i][j]);
          if (Length_ref > 0.0) Bond_ff[i][j]/=Length_ref;
        }
      }
    }
    else  {
      for (i=0; i<Ncomp; i++){
         if (Mix_type==0) {jmin=i; jmax=i+1;}
         else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
         for(j=jmin;j<jmax;j++) Bond_ff[i][j] = 0.0;
      }
      fprintf(fp2,"n/a: bonds only needed for polymer runs");
    }
  }
  MPI_Bcast(Bond_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

                       /* WALL-WALL PARAMS */

                       /* Density of atoms in the surfaces */
  if (Proc==0) {
    read_junk(fp,fp2);
      for (i=0; i<Nwall_type; i++){
        fscanf(fp,"%lf", &Rho_w[i]);
        fprintf(fp2,"%f  ",Rho_w[i]);
        if (Density_ref > 0.0) Rho_w[i] /=Density_ref;
      }
  }
  MPI_Bcast(Rho_w,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
                       /*Sigma_w or Sigma_ww*/
  if (Proc==0) {
    read_junk(fp,fp2);
      for (i=0; i<Nwall_type; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Nwall_type;}
        for (j=jmin; j<jmax; j++){
	  if (Mix_type==1){ fscanf(fp,"%lf",&Sigma_ww[i][j]);
	                    fprintf(fp2,"%f  ",Sigma_ww[i][j]);
                            if (Length_ref > 0.0) Sigma_ww[i][j]/=Length_ref;
                            if (Surface_type[i] == atomic_centers && j==i) 
                	         WallParam[i] = Sigma_ww[i][i]/2.0;
                          }
          else            { fscanf(fp,"%lf",&Sigma_w[i]);
	                    fprintf(fp2,"%f  ",Sigma_w[i]);
                            if (Length_ref > 0.0) Sigma_w[i]/=Length_ref;
                            if (Surface_type[i] == atomic_centers) 
                	         WallParam[i] = Sigma_w[i]/2.0;
                          }
        }
      }
  }
    if (Mix_type==1) MPI_Bcast(Sigma_ww,NWALL_MAX_TYPE*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    else MPI_Bcast(Sigma_w,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

                       /*Eps_w or Eps_ww*/
  if (Proc==0) {
    read_junk(fp,fp2);
      for (i=0; i<Nwall_type; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Nwall_type;}
        for (j=jmin; j<jmax; j++){
	  if (Mix_type==1){ fscanf(fp,"%lf",&Eps_ww[i][j]);
	                    fprintf(fp2,"%f  ",Eps_ww[i][j]);
                            if (Temp > 0.0) Eps_ww[i][j]/=Temp;
                          }
          else            { fscanf(fp,"%lf",&Eps_w[i]);
	                    fprintf(fp2,"%f  ",Eps_w[i]);
                            if (Temp > 0.0) Eps_w[i]/=Temp;
                          }
        }
      }
  }
    if (Mix_type==1) MPI_Bcast(Eps_ww,NWALL_MAX_TYPE*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    else MPI_Bcast(Eps_w,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

                       /*Cut_ww */
  if (Proc==0) {
    read_junk(fp,fp2);
      for (i=0; i<Nwall_type; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Nwall_type;}
        for (j=jmin; j<jmax; j++){
	   fscanf(fp,"%lf",&Cut_ww[i][j]);
	                    fprintf(fp2,"%f  ",Cut_ww[i][j]);
                            if (Length_ref > 0.0) Cut_ww[i][j]/=Length_ref;
        }
      }
  }
    MPI_Bcast(Cut_ww,NWALL_MAX_TYPE*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);


                       /* WALL-FLUID PARAMS -- ONLY IF MIX_TYPE == 1 */

  if (Mix_type == 1){

     if (Proc==0) {
       read_junk(fp,fp2);
                       /* Sigma_wf */
         for (i=0; i<Ncomp; i++){
           for (j=0; j<Nwall_type; j++){
   	  fscanf(fp,"%lf",&Sigma_wf[i][j]);
   	  fprintf(fp2,"%f  ",Sigma_wf[i][j]);
             if (Length_ref > 0.0) Sigma_wf[i][j]/=Length_ref;
           }
         }
                       /* Eps_wf */
       read_junk(fp,fp2);
         for (i=0; i<Ncomp; i++){
           for (j=0; j<Nwall_type; j++){
   	  fscanf(fp,"%lf",&Eps_wf[i][j]);
   	  fprintf(fp2,"%f  ",Eps_wf[i][j]);
             if (Temp > 0.0) Eps_wf[i][j]/=Temp;
           }
         }
                       /* Cut_wf */
       read_junk(fp,fp2);
         for (i=0; i<Ncomp; i++){
           for (j=0; j<Nwall_type; j++){
   	  fscanf(fp,"%lf",&Cut_wf[i][j]);
   	  fprintf(fp2,"%f  ",Cut_wf[i][j]);
             if (Length_ref > 0.0) Cut_wf[i][j]/=Length_ref;
           }
         }
     }
    MPI_Bcast(Sigma_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Eps_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Cut_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  else{
    if (Proc==0){
         read_junk(fp,fp2);
         fprintf(fp2,"\n  USING L_B Mixing Rules --- WALL-FLUID INTERACTIONS COMPUTED BY CODE\n");
         fprintf(fp2,"........MANUAL INPUT DOES NOT APPLY \n");
         fprintf(fp2,"skipping this parameter  ");
         for (i=0; i<2; i++){
                read_junk(fp,fp2);
                fprintf(fp2,"skipping this parameter  ");
         }
    }
  }

  /* Read and set up Polymer parameters - only read if stencil is turned on */
  Nbonds=0;
  for (icomp=0;icomp<Ncomp;icomp++){
       Nseg_type[icomp]=1;  /* if no polymer we just set to 1 */
       Unk2Comp[icomp]=icomp;
  }

  if (Type_poly != NONE){

    if (Proc==0) {
      printf("\n");
      read_junk(fp,fp2);
      fscanf(fp,"%d",&Npol_comp);
      fprintf(fp2,"%d",Npol_comp);
    }
    MPI_Bcast(&Npol_comp,1,MPI_INT,0,MPI_COMM_WORLD);
    if (Proc==0) {
      read_junk(fp,fp2);
      for (i=0; i<Npol_comp; ++i){
	fscanf(fp,"%d",&Nblock[i]);
	fprintf(fp2,"%d",Nblock[i]);
	if (Nblock[i] > NBLOCK_MAX) {
	  if (Proc==0) printf("Must increase NBLOCK_MAX");
	  Nblock[i] = NBLOCK_MAX;
	}
      }
    }
    MPI_Bcast(Nblock,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
    if(Proc==0) {
      read_junk(fp,fp2);
      for (pol_number=0; pol_number<Npol_comp; ++pol_number){
	Nmer[pol_number] = 0;
	for (i=0; i<Nblock[pol_number]; ++i){
	  fscanf(fp,"%d", &block[pol_number][i]);
	  fprintf(fp2,"%d  ",block[pol_number][i]);
	  Nmer[pol_number] += block[pol_number][i];
	}
	printf("Length (Nmer) of polymer %d = %d\n",pol_number,
	       Nmer[pol_number]);
      }
    }
    MPI_Bcast(Nmer,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
    if (Proc==0) {
      Ntype_mer = 0;
      read_junk(fp,fp2);
      seg_tot=0;
      for (i=0; i<NBLOCK_MAX; ++i) Nmer_t_total[i]=0;
      for (pol_number=0; pol_number<Npol_comp; ++pol_number){
	seg = 0;
	for (i=0; i<NBLOCK_MAX; ++i)  Nmer_t[pol_number][i] = 0; 
	for (i=0; i<Nblock[pol_number]; ++i){
	  fscanf(fp,"%d", &block_type[i]);
	  fprintf(fp2,"%d  ",block_type[i]);
	  Nmer_t[pol_number][block_type[i]] += block[pol_number][i];
	  Nmer_t_total[block_type[i]] += block[pol_number][i];
	  for (j=0; j<block[pol_number][i]; j++) {
	    Type_mer[pol_number][seg] = block_type[i];
            SegChain2SegAll[pol_number][seg]=seg_tot;
            seg++; seg_tot++;
          }
	  if (block_type[i] > Ntype_mer) Ntype_mer = block_type[i];
	}
      }
      Ntype_mer++;
      printf("Number of segment types (Ntype_mer) = %d\n",Ntype_mer);
    }
    MPI_Bcast(&Ntype_mer,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Nmer_t,NCOMP_MAX*NBLOCK_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Nmer_t_total,NBLOCK_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(SegChain2SegAll,NCOMP_MAX*NMER_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Type_mer,NCOMP_MAX*NMER_MAX,MPI_INT,0,MPI_COMM_WORLD);


    if (Proc==0) {
      read_junk(fp,fp2);
      fscanf(fp,"%s", poly_file);
      fprintf(fp2,"%s  ",poly_file);
      if( (fp4  = fopen(poly_file,"r")) == NULL) {
	printf("Can't open file %s\n", poly_file);
	exit(1);
      }
    }

    /* now read in the bonds linking the various segments - note that there
       may be more than two bonds per site */
  
    nseg=nmer_max=0;
    for (pol_number=0; pol_number<Npol_comp; pol_number++){
         nseg += Nmer[pol_number];
         if (Nmer[pol_number] > nmer_max) nmer_max=Nmer[pol_number];
    }
    Unk_to_Poly = (int *) array_alloc (1, NBOND_MAX*nseg, sizeof(int));
    Unk_to_Seg = (int *) array_alloc (1, NBOND_MAX*nseg, sizeof(int));
    Unk_to_Bond = (int *) array_alloc (1, NBOND_MAX*nseg, sizeof(int));

    Nbond = (int **) array_alloc (2, Npol_comp,nmer_max,sizeof(int));
    Nbonds_SegAll = (int *) array_alloc (1, NMER_MAX,sizeof(int));
    Bonds = (int ***) array_alloc (3, Npol_comp,nmer_max,NBOND_MAX,sizeof(int));
    Bonds_SegAll = (int **) array_alloc (2, NMER_MAX,NBOND_MAX,sizeof(int));
    pol_sym_tmp = (int ***) array_alloc (3, Npol_comp,nmer_max,NBOND_MAX,sizeof(int));
    Poly_to_Unk = (int ***) array_alloc (3, Npol_comp,nmer_max,NBOND_MAX,sizeof(int));
    Poly_to_Unk_SegAll = (int **) array_alloc (2, NMER_MAX,NBOND_MAX,sizeof(int));
    Pol_Sym = (int *) array_alloc (1, nseg*NBOND_MAX,sizeof(int));
    BondAll_to_isegAll = (int *) array_alloc (1, nseg*NBOND_MAX,sizeof(int));
    BondAll_to_ibond = (int *) array_alloc (1, nseg*NBOND_MAX,sizeof(int));
    nbond_tot = (int *) array_alloc (1, Npol_comp, sizeof(int));

    nbond_all = 0; 
    Nbonds=0;
    Nseg_tot=0;
    seg_tot=0;
    for (pol_number=0; pol_number<Npol_comp; ++pol_number){
      Nseg_tot += Nmer[pol_number];
      nbond_tot[pol_number]=0;
      nunk=0; 
      for (iseg=0; iseg<Nmer[pol_number]; iseg++){
        end_count=0;
	if (Proc==0) fscanf(fp4,"%d", &Nbond[pol_number][iseg]);
	MPI_Bcast(&Nbond[pol_number][iseg],1,MPI_INT,0,MPI_COMM_WORLD);
        Nbonds_SegAll[seg_tot]=0;

	for (ibond=0; ibond<Nbond[pol_number][iseg]; ibond++){
	  if (Proc==0) {
	    fscanf(fp4,"%d  %d", &Bonds[pol_number][iseg][ibond],&pol_sym_tmp[pol_number][iseg][ibond]);
	    fprintf(fp2,"%d  ",Bonds[pol_number][iseg][ibond]);
	  }
	  MPI_Bcast(&Bonds[pol_number][iseg][ibond],1,MPI_INT,0,MPI_COMM_WORLD);
	  MPI_Bcast(&pol_sym_tmp[pol_number][iseg][ibond],1,MPI_INT,0,MPI_COMM_WORLD);
          if (Type_poly!=WTC || (Type_poly==WTC && Bonds[pol_number][iseg][ibond] != -1)){
                                  /* note we don't want to include ends for the WTC polymers!*/
	    Unk_to_Poly[nbond_all] = pol_number;
  	    Unk_to_Seg[nbond_all]  = iseg;
	    Unk_to_Bond[nbond_all] = ibond;
	    Poly_to_Unk[pol_number][iseg][ibond] = nunk;
            Bonds_SegAll[seg_tot][Nbonds_SegAll[seg_tot]]=Bonds[pol_number][iseg][ibond]+SegChain2SegAll[pol_number][0];
	    Poly_to_Unk_SegAll[seg_tot][Nbonds_SegAll[seg_tot]] = nunk;
	    Pol_Sym[nbond_all]=pol_sym_tmp[pol_number][iseg][ibond];
            BondAll_to_isegAll[nbond_all]=seg_tot;
            BondAll_to_ibond[nbond_all]=ibond;
	    nbond_all++;
	    nunk++;
            Nbonds++; 
            Nbonds_SegAll[seg_tot]++;
          }
          else if (Type_poly==WTC && Bonds[pol_number][iseg][ibond] == -1){
              end_count++;
          }
	}
	nbond_tot[pol_number] += (Nbond[pol_number][iseg]-end_count);
        Unk2Comp[seg_tot]=Type_mer[pol_number][iseg];
        seg_tot++;
      }
    }
    for (icomp=0;icomp<Ncomp;icomp++) Nseg_type[icomp]=0;
    for (iseg=0;iseg<Nseg_tot;iseg++) Nseg_type[Unk2Comp[iseg]]++;
    if (Proc==0){
       fprintf(fp2,"\n********************\n BOND DETAILS \n **********************\n");
       fprintf(fp2,"\t total number of bonds is %d\n",Nbonds);
       for (ibond=0;ibond<Nbonds; ibond++){
           fprintf(fp2,"Unk_to_Poly[ibond=%d]=%d Unk_to_Seg[]=%d Unk_to_Bond[]=%d\n",
            ibond,Unk_to_Poly[ibond],Unk_to_Seg[ibond],Unk_to_Bond[ibond]);
       }
       for (pol_number=0; pol_number<Npol_comp; ++pol_number){
          for (iseg=0;iseg<Nmer[pol_number];iseg++){
              for (ibond=0;ibond<Nbond[pol_number][iseg];ibond++){
                  fprintf(fp2,"Poly_to_Unk[%d][%d][%d]=%d\n",
                     pol_number, iseg,ibond,Poly_to_Unk[pol_number][iseg][ibond]);
              }
          }
       }
       for (pol_number=0; pol_number<Npol_comp; ++pol_number){
          for (iseg=0;iseg<Nmer[pol_number];iseg++){
              printf("SegChain2SegAll[%d][%d]=%d\n",
                   pol_number, iseg,SegChain2SegAll[pol_number][iseg]);
          }
       }
       printf("Total Number of segments in the problem=%d\n",Nseg_tot);
       for (iseg=0;iseg<Nseg_tot;iseg++){
           printf("Nbonds_SegAll[%d]=%d",iseg,Nbonds_SegAll[iseg]);
           printf("\t seg %d is bonded to...",iseg);
           for(ibond=0;ibond<Nbonds_SegAll[iseg];ibond++)
                printf("%d  ",Bonds_SegAll[iseg][ibond]); 
           printf("\n");
       }
       fprintf(fp2,"****************\n END BOND DETAILS \n **********************\n");
    }
    
    if (Proc==0) fclose(fp4);

    if (Sten_Type[POLYMER_CR]){  /*POLYMER INPUT FOR ONLY CMS FUNCTIONAL */
    /* set start value of Geqns for each of the polymers in the system.  
       It is necessary to account for Ncomp Boltz and Ncomp Rho eqns */
   
       Ngeqn_tot=0; 
       for (pol_number=0; pol_number<Npol_comp; ++pol_number){
          Geqn_start[pol_number] = 0;
          Ngeqn_tot += (nbond_tot[pol_number]);
          for (pol_num2=0; pol_num2<pol_number; pol_num2++)
              Geqn_start[pol_number] += (nbond_tot[pol_num2]);
       }
       safe_free((void *)  &nbond_tot); 
       if (Proc==0) printf("The total number of g equations will be %d\n",Ngeqn_tot);
       for (pol_number=0; pol_number<Npol_comp; ++pol_number)
       if (Proc==0) printf("The start unknown for polymer %d is %d \n",
                                    pol_number,Geqn_start[pol_number]);
    
       if (Proc==0) {
          read_junk(fp,fp2);
          fscanf(fp,"%d",&Ncr_files);
          fscanf(fp,"%lf", &Crfac);
          fscanf(fp,"%s", Cr_file);
          fprintf(fp2,"%d  %9.6f  %s ", Ncr_files,Crfac,Cr_file);
          if (Ncr_files>=2){ 
                fscanf(fp,"%s",Cr_file2);
                fprintf(fp2,"  %s ",Cr_file2);
          }
          if (Ncr_files>=3){ 
                fscanf(fp,"%s",Cr_file3);
                fprintf(fp2,"  %s ",Cr_file3);
          }
          if (Ncr_files>=4){ 
                fscanf(fp,"%s",Cr_file4);
                fprintf(fp2,"  %s ",Cr_file4);
          }
          read_junk(fp,fp2);
/*        if ( fabs(Crfac+1.0)<1.e-8){*/
            for (i=0;i<Ncr_files-2;i++){
               fscanf(fp,"%lf", &Cr_break[i]);
               fprintf(fp2,"%f  ",Cr_break[i]);
            }
        /*}
          else fprintf(fp2,"n/a: not doing automated linear interpolation");*/
          read_junk(fp,fp2);
          for  (icomp=0; icomp<Ncomp; icomp++){
             for  (jcomp=0; jcomp<Ncomp; jcomp++){
                fscanf(fp,"%lf", &Cr_rad_hs[icomp][jcomp]);
                fprintf(fp2,"%f  ",Cr_rad_hs[icomp][jcomp]);
             }
          }
       }
       MPI_Bcast(Cr_rad_hs,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
       MPI_Bcast(Cr_break,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
       MPI_Bcast(&Ncr_files,1,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&Crfac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    /* Note: the value of Cr_rad_hs may get reset to Cutoff_ff if Ipot_ff_n=2 
       see setup_polymer_cr in dft_main.c */
       if (Proc==0) {
          read_junk(fp,fp2);
          fscanf(fp,"%lf", &Gauss_a);
          if (Proc==0) fprintf(fp2,"%f  ",Gauss_a);
       }
       MPI_Bcast(&Gauss_a,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
       r = 3./(2.*PI*Gauss_a*Gauss_a);
       Gauss_k = pow(r,1.5);

/* This is a debugging tool used to enter a real and an
   integer to try various things out.  Generally, we would
   say at iter # Bupdate_iters set a certain number to 
   Bupdate_fact.  Since this is useless for the average
   user, these numbers are just set to nonsense here in
   the input file, if they are needed, just uncomment these
   lines, and add the input to the input file !
*/
/*    
       if (Proc==0) {
         read_junk(fp,fp2);
         fscanf(fp,"%lf", &Bupdate_fact);
         fprintf(fp2,"%f  ",Bupdate_fact);
         fscanf(fp,"%d", &Bupdate_iters);
         fprintf(fp2,"%d  ",Bupdate_iters);
       }
       MPI_Bcast(&Bupdate_fact,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
       MPI_Bcast(&Bupdate_iters,1,MPI_INT,0,MPI_COMM_WORLD);
*/
       Bupdate_fact = 0.0;
       Bupdate_iters = -1;
    } /* ends the section on Geqns that is only used if POLYMER_CR stencil is on */
    else{
       if (Proc==0) {
          read_junk(fp,fp2);
          fprintf(fp2,"\n NO LIQUID STATE INPUT FOR WERTHEIM-TRIPATHI-CHAPMAN RUN\n");
          fprintf(fp2,"not read   ");
          for (i=0; i<3; i++) {
             read_junk(fp,fp2);
             fprintf(fp2,"not read   ");
          }       
       }       
    }
  }     /* POLYMER INPUT FOR ALL KINDS OF POLYMERS */
  else{
    if (Proc==0) {
      read_junk(fp,fp2);
      fprintf(fp2,"\n POLYMER INPUT NOT RELEVENT FOR THIS RUN\n");
      fprintf(fp2,"not read   ");
      for (i=0; i<8; i++) {
             read_junk(fp,fp2);
             fprintf(fp2,"not read   ");
      }
    }
  }
  /* end of polymer input */

  /* Read in  Semi-Permeable Surface parameters */
  if (Nwall_type >0) Lsemiperm = (int **) array_alloc (2, Nwall_type,Ncomp,sizeof(int));

  lzeros=FALSE;
  ltrues=FALSE;
  if (Proc==0) {
    read_junk(fp,fp2);
    for (i=0; i<Nwall_type; i++){
      for(j=0; j<Ncomp; j++){
        if (!lzeros  && ! ltrues){
           fscanf(fp,"%d ", &Lsemiperm[i][j]);
           fprintf(fp2,"%d  ",Lsemiperm[i][j]);
        }
        if (i==0 && j==0){
           if (Lsemiperm[i][j]==-2) lzeros=TRUE;
           else if (Lsemiperm[i][j]==-1) ltrues=TRUE;
        }
        if (lzeros) Lsemiperm[i][j]=0;
        else if (ltrues) Lsemiperm[i][j]=1;
      }
    }
  }
  if (Nwall_type >0)
  MPI_Bcast(*Lsemiperm,Nwall_type*Ncomp,MPI_INT,0,MPI_COMM_WORLD);


  if (Nwall_type>0)
  Vext_membrane = (double **) array_alloc (2, Nwall_type,Ncomp,sizeof(double));
  if (Proc==0) {
    read_junk(fp,fp2);
    for (i=0; i<Nwall_type; i++){
      for(j=0; j<Ncomp; j++){
        if (lzeros){
           Vext_membrane[i][j]=0.0;
        }
        else{
          fscanf(fp,"%lf", &Vext_membrane[i][j]);
          fprintf(fp2,"%f  ",Vext_membrane[i][j]);
        }
        if (Temp > 0.0) Vext_membrane[i][j] /= Temp;
      }
    }
    for (i=0; i<Nwall_type; i++){
      for(j=0; j<Ncomp; j++)
	if (! Lsemiperm[i][j]) Vext_membrane[i][j] = VEXT_MAX;
    }
  }
  if (Nwall_type > 0)
  MPI_Bcast(*Vext_membrane,Nwall_type*Ncomp,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /* State Point Parameters */


  if (Proc==0) {
    read_junk(fp,fp2);
    if (Type_poly==NONE || Ntype_mer == 1){
      for (icomp=0; icomp<Ncomp; ++icomp){
        fscanf(fp,"%lf", &Rho_b[icomp]);
        fprintf(fp2,"%f  ",Rho_b[icomp]);
        if (Density_ref > 0.0) Rho_b[icomp] /= Density_ref;
      }
    }
    else{
      /*first scan in polymer densities */
      for (i=0; i<Npol_comp; i++){
          fscanf(fp,"%lf", &rho_tmp[i]);
          if (Density_ref > 0.0) rho_tmp[i] /= Density_ref;
      }
      /* now scan in solvent - atom densities */
      for (icomp=Ntype_mer; icomp<Ncomp; ++icomp){
        fscanf(fp,"%lf", &Rho_b[icomp]);
        if (Density_ref > 0.0) Rho_b[icomp] /= Density_ref;
      }
      for (i=0; i<Ntype_mer; i++){
        Rho_b[i] = 0.;
        for (j=0; j<Npol_comp; j++) {
          Rho_b[i] += (double)Nmer_t[j][i]*rho_tmp[j]/(double)Nmer[j];
        }
        fprintf(fp2,"%f  ",Rho_b[i]);
      }
      for (icomp=Ntype_mer; icomp<Ncomp; ++icomp)
        fprintf(fp2,"%f  ",Rho_b[icomp]);
    }
  }
  MPI_Bcast(Rho_b,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /* Read in Charged surface parameters */
  Ipot_wf_c = 0;
  if (Type_coul != -1) {
    if (Proc==0) {
      read_junk(fp,fp2);
      for (iwall_type=0; iwall_type <  Nwall_type; ++iwall_type){
	fscanf(fp,"%d", &Type_bc_elec[iwall_type]);
	fprintf(fp2,"%d  ",Type_bc_elec[iwall_type]);
        if (Type_bc_elec[iwall_type] != 0) Ipot_wf_c = COULOMB;
      }
    }
    if (Nwall_type>0) MPI_Bcast(Type_bc_elec,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);

    /* Make unit corrections to the Elec_param_w array if needed */
    for (iwall=0; iwall <  Nwall; ++iwall){
        if(Type_bc_elec[WallType[iwall]]==1 && Potential_ref>0.0) 
                                 Elec_param_w[iwall]/=Potential_ref;
        if(Type_bc_elec[WallType[iwall]]==2 && Length_ref>0.0) 
                                 Elec_param_w[iwall]*=(Length_ref*Length_ref);
    }

    if (Proc==0) {
      read_junk(fp,fp2);
      fscanf(fp,"%d", &Nlocal_charge);
      fprintf(fp2,"%d  ",Nlocal_charge);
      if (Nlocal_charge > 0) Ipot_wf_c=COULOMB;
    }
    MPI_Bcast(&Nlocal_charge,1,MPI_INT,0,MPI_COMM_WORLD);
    ncharge=0;
    if      (Nlocal_charge > 0)  ncharge = Nlocal_charge;
    else if (Nlocal_charge < 0)  ncharge = 2;
    
    if (Proc==0) read_junk(fp,fp2);
    if (Nlocal_charge == 0)  {
      if (Proc==0) {
	fprintf(fp2,"n/a: Nlocal_charge=0");
	read_junk(fp,fp2);
	fprintf(fp2,"n/a: Nlocal_charge=0");
	read_junk(fp,fp2);
	fprintf(fp2,"n/a: Nlocal_charge=0");
      }
    }
    else{
      Charge      = (double *) array_alloc (1, ncharge, sizeof(double));
      Charge_Diam = (double *) array_alloc (1, ncharge, sizeof(double));
      Charge_x = (double **) array_alloc (2, Ndim,ncharge,sizeof(double));
      if (Proc==0) {
	for(i=0; i<ncharge; i++){
	  fscanf(fp,"%lf ", &Charge[i]);
	  if (Proc==0) fprintf(fp2,"%f  ",Charge[i]);
	}
      }
      MPI_Bcast(Charge,ncharge,MPI_DOUBLE,0,MPI_COMM_WORLD);
      if (Proc==0) {
	read_junk(fp,fp2);
	for(i=0; i<ncharge; i++){
	  fscanf(fp,"%lf ", &Charge_Diam[i]);
	  fprintf(fp2,"%f  ",Charge_Diam[i]);
	}
      }
      MPI_Bcast(Charge_Diam,ncharge,MPI_DOUBLE,0,MPI_COMM_WORLD);
      if (Proc==0) {
	read_junk(fp,fp2);
	for(i=0; i<ncharge; i++){
	  for (idim=0; idim <  Ndim; ++idim){
	    fscanf(fp,"%lf", &Charge_x[idim][i]);
	    fprintf(fp2,"%f  ",Charge_x[idim][i]);
	  }
	}
      }
      MPI_Bcast(*Charge_x,Ndim*ncharge,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    logical=FALSE;
    for (iwall_type=0; iwall_type<Nwall_type; iwall_type++)
      if (Type_bc_elec[iwall_type] == ATOMIC_CHARGE) logical=TRUE;

    if (Proc==0)  read_junk(fp,fp2);
    if (ncharge !=0 || logical){
      if (Proc==0) {
	fscanf(fp,"%d  %d", &Charge_type_atoms, &Charge_type_local);
	fprintf(fp2,"%d  %d ",Charge_type_atoms,Charge_type_local);
      }
      MPI_Bcast(&Charge_type_atoms,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&Charge_type_local,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    else if (Proc==0)  fprintf(fp2,"n/a:  ncharge=0 and walls are not atoms");
  }
  else{
    if (Proc==0) {
      read_junk(fp,fp2);
      fprintf(fp2,"CHARGED SURFACE INPUT NOT RELEVENT FOR THIS RUN\n");
      fprintf(fp2,"n/a");
      for (i=0; i<5; i++){
	read_junk(fp,fp2);
	fprintf(fp2,"n/a");
      }
    }
  }
  MPI_Bcast(&Ipot_wf_c,1,MPI_INT,0,MPI_COMM_WORLD);
  
   /* Dielectric constant choices */

  if (Ipot_ff_c == COULOMB ) {
    if (Proc==0) {
      read_junk(fp,fp2);
      fscanf(fp,"%d  ",&Type_dielec);
      fprintf(fp2,"%d  ",Type_dielec);
    }
    MPI_Bcast(&Type_dielec,1,MPI_INT,0,MPI_COMM_WORLD);
    if (Length_ref > 0 && Dielec_ref > 0 && Temp > 0){
       Temp_elec = 4.0*PI*KBOLTZ*Temp*Dielec_ref*EPSILON_0*Length_ref*1.e-10/(E_CONST*E_CONST);
       if (Proc==0) printf("\t WARNING :: CHECK THAT Length_ref WAS GIVEN IN ANGSTROM UNITS ... CHARGED SYSTEM\n");
    }
    else /* make some assumptions */
       Temp_elec = 4*PI*KBOLTZ*298.0*KAPPA_H2O*EPSILON_0*4.25e-10/(E_CONST*E_CONST);

    if (Proc==0) printf("\t Temp_elec=%9.6f\n",Temp_elec);

    if (Proc==0) read_junk(fp,fp2);
    if (Type_dielec != DIELEC_WF_PORE){
      if (Proc==0) {
	fscanf(fp,"%lf",&Dielec_bulk);
	fprintf(fp2,"%f  ",Dielec_bulk);
        if (Dielec_ref > 0.0) Dielec_bulk /= Dielec_ref;
      }
      MPI_Bcast(&Dielec_bulk,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    else {
      if (Proc==0) {
	fscanf(fp,"%lf %lf %lf",&Dielec_bulk, &Dielec_pore, &Dielec_X);
	fprintf(fp2,"%f  %f  %f ",Dielec_bulk,Dielec_pore,Dielec_X);
        if (Dielec_ref > 0.0) Dielec_bulk /= Dielec_ref;
        if (Dielec_ref > 0.0) Dielec_pore /= Dielec_ref;
        if (Length_ref > 0.0) Dielec_X /= Length_ref;
      }
      MPI_Bcast(&Dielec_bulk,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&Dielec_pore,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&Dielec_X,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    
    if (Proc==0) read_junk(fp,fp2);
    if (Nwall_type > 0) {
      Dielec_wall = (double *) array_alloc (1, Nwall_type,sizeof(double));
      if (Proc==0) {
	for (iwall_type=0; iwall_type<Nwall_type; iwall_type++) {
	  fscanf(fp,"%lf", &Dielec_wall[iwall_type]);
	  fprintf(fp2,"%f  ",Dielec_wall[iwall_type]);
          if (Dielec_ref > 0.0) Dielec_wall[iwall_type]/=Dielec_ref;
	}
      }
      MPI_Bcast(Dielec_wall,Nwall_type,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    else if (Proc==0)  fprintf(fp2,"n/a");
  }
  else{
    if (Proc==0) {
      read_junk(fp,fp2);
      fprintf(fp2,"CHARGED FLUID (DIELECTRIC) INPUT NOT RELEVENT FOR THIS RUN\n");
      fprintf(fp2,"n/a");
      for (i=0; i<2; i++){
	read_junk(fp,fp2);
	fprintf(fp2,"n/a");
      }
    }
  }
  
  /* READ IN STEADY STATE BOUNDARY CONDITIONS */
  L1D_bc = FALSE;
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d  ", &Lsteady_state);
    fprintf(fp2,"%d  ",Lsteady_state);
    Linear_transport=0;  
                        /* only want generalized Fick's Law */
                        /* can do J=-Dgrad mu, but doesn't work well */
  }
  MPI_Bcast(&Lsteady_state,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Linear_transport,1,MPI_INT,0,MPI_COMM_WORLD);

  /* ALF: temporary fix for polymer initial guesses */
  if (Lsteady_state || Sten_Type[POLYMER_CR]) {
    if (Proc==0) {
      read_junk(fp,fp2);
      fscanf(fp,"%d  %d  %lf", &Grad_dim, &L1D_bc,&X_1D_bc );
      fprintf(fp2,"%d  %d  %f  ",Grad_dim, L1D_bc,X_1D_bc);
      if (Length_ref > 0.0) X_1D_bc /= Length_ref;
    }
    MPI_Bcast(&Grad_dim,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&L1D_bc,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&X_1D_bc,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (Proc==0) {
      read_junk(fp,fp2);
      fscanf(fp,"%lf", &X_const_mu);
      fprintf(fp2,"%f  ",X_const_mu);
      if (Length_ref > 0.0) X_const_mu /= Length_ref;
    }
    MPI_Bcast(&X_const_mu,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (Proc==0) {
      read_junk(fp,fp2);
      fscanf(fp,"%d  %d", &Geom_flag, &Nseg_IC);
      fprintf(fp2,"%d  %d  ",Geom_flag,Nseg_IC);
    }
    MPI_Bcast(&Geom_flag,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&Nseg_IC,1,MPI_INT,0,MPI_COMM_WORLD);
    Pore_rad_L_IC = (double *) array_alloc (1, Nseg_IC,sizeof(double));
    Pore_rad_R_IC = (double *) array_alloc (1, Nseg_IC,sizeof(double));
    Lseg_IC       = (double *) array_alloc (1, Nseg_IC,sizeof(double));
    if (Proc==0) {
      read_junk(fp,fp2);
      for (i=0; i<Nseg_IC; i++){
	fscanf(fp,"%lf  %lf  %lf", &Pore_rad_L_IC[i], &Pore_rad_R_IC[i], &Lseg_IC[i]);
	fprintf(fp2,"%f  %f  %f  ", Pore_rad_L_IC[i],Pore_rad_R_IC[i],Lseg_IC[i]);
        if (Length_ref > 0.0){
            Pore_rad_L_IC[i]/=Length_ref;
            Pore_rad_R_IC[i]/=Length_ref;
            Lseg_IC[i]/=Length_ref;
        }
      }
    }
    MPI_Bcast(Pore_rad_L_IC,Nseg_IC,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Pore_rad_R_IC,Nseg_IC,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Lseg_IC,Nseg_IC,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (Proc==0) {
      read_junk(fp,fp2);
      for (icomp=0; icomp<Ncomp; ++icomp){
	fscanf(fp,"%lf", &Rho_b_LBB[icomp]);
	fprintf(fp2,"%f  ",Rho_b_LBB[icomp]);
        if (Density_ref > 0.0) Rho_b_LBB[icomp] /= Density_ref;
      }
    }
    MPI_Bcast(Rho_b_LBB,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (Proc==0) {
      read_junk(fp,fp2);
      for (icomp=0; icomp<Ncomp; ++icomp){
	fscanf(fp,"%lf", &Rho_b_RTF[icomp]);
	fprintf(fp2,"%f  ",Rho_b_RTF[icomp]);
        if (Density_ref > 0.0) Rho_b_RTF[icomp] /= Density_ref;
      }
    }
    MPI_Bcast(Rho_b_RTF,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (Proc==0) {
      read_junk(fp,fp2);
      for (icomp=0; icomp<Ncomp; ++icomp){
	fscanf(fp,"%lf", &D_coef[icomp]);
	fprintf(fp2,"%f  ",D_coef[icomp]);
      }
    }
    MPI_Bcast(D_coef,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if (Ipot_ff_c == COULOMB) {
      if (Proc==0) {
	read_junk(fp,fp2);
	fscanf(fp,"%lf", &Elec_pot_LBB);
	fprintf(fp2,"%f  ",Elec_pot_LBB);
        if (Potential_ref > 0.0) Elec_pot_LBB /= Potential_ref;
      }
      MPI_Bcast(&Elec_pot_LBB,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      if (Proc==0) {
	fscanf(fp,"%lf", &Elec_pot_RTF);
	fprintf(fp2,"%f  ",Elec_pot_RTF);
        if (Potential_ref > 0.0) Elec_pot_RTF /= Potential_ref;
      }
      MPI_Bcast(&Elec_pot_RTF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    else{
      if (Proc==0) {
	read_junk(fp,fp2);
	fprintf(fp2,"n/a");
      }
    }
    if (Proc==0) {
      read_junk(fp,fp2);
      fscanf(fp,"%lf ", &Velocity);
      fprintf(fp2,"%f  ",Velocity);
    }
    MPI_Bcast(&Velocity,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  else{
    if (Proc==0) {
      read_junk(fp,fp2);
      fprintf(fp2,"TRANSPORT INPUT NOT RELEVENT FOR THIS RUN\n");
      fprintf(fp2,"n/a");
      for (i=0; i<8; i++) {
	read_junk(fp,fp2);
	fprintf(fp2,"n/a");
      }
    }
  }

  /* ***********end of physics input .... now run control and numerics input ********/ 

  /* Run control parameters: initial guess type etc.*/

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Iliq_vap);
    fprintf(fp2,"%d",Iliq_vap);
  }
  MPI_Bcast(&Iliq_vap,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Iguess1);
    fprintf(fp2,"%d   ",Iguess1);
  }
  MPI_Bcast(&Iguess1,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Iguess1==STEP_PROFILE || (Iguess1>=CHOP_RHO && Iguess1<= CHOP_RHO_STEP)){
    if (Proc==0) {
      read_junk(fp,fp2);
      if (Iguess1 != STEP_PROFILE) Nsteps=1;
      else   fscanf(fp,"%d",&Nsteps);
      fprintf(fp2,"%d   ",Nsteps);
    }
    MPI_Bcast(&Nsteps,1,MPI_INT,0,MPI_COMM_WORLD);
    if (Proc==0) {
      read_junk(fp,fp2);
      for (i=0; i<Nsteps; ++i){
	fscanf(fp,"%d", &Orientation_step[i]);
	fprintf(fp2,"%d  ",Orientation_step[i]);
      }
    }
    MPI_Bcast(Orientation_step,NSTEPS_MAX,MPI_INT,0,MPI_COMM_WORLD);
    if (Proc==0) {
      read_junk(fp,fp2);
      for (i=0; i<Nsteps; ++i){
	fscanf(fp,"%lf", &Xstart_step[i]);
	fprintf(fp2,"%f  ",Xstart_step[i]);
        Xstart_step[i]-=1.0e-4; /* prevent roundoff errors later*/
      }
    }
    MPI_Bcast(Xstart_step,NSTEPS_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (Proc==0) {
      read_junk(fp,fp2);
      for (i=0; i<Nsteps; ++i){
	fscanf(fp,"%lf", &Xend_step[i]);
	fprintf(fp2,"%f  ",Xend_step[i]);
        Xend_step[i]+=1.e-4; /* prevent roundoff errors later*/
      }
    }
    MPI_Bcast(Xend_step,NSTEPS_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (Proc==0) {
      read_junk(fp,fp2);
      for (icomp=0; icomp<Ncomp; ++icomp){
          for (i=0; i<Nsteps; ++i){
   	    fscanf(fp,"%lf", &Rho_step[icomp][i]);
	    fprintf(fp2,"%f  ",Rho_step[icomp][i]);
          }
      }
    }
    MPI_Bcast(Rho_step,NSTEPS_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  else{
    if (Proc==0) for (i=0;i<5;i++) read_junk(fp,fp2);
  }

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Restart);
    fprintf(fp2,"%d",Restart);
  }
  MPI_Bcast(&Restart,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
	fscanf(fp,"%lf",&Rho_max);
	fprintf(fp2,"%f  ",Rho_max);
  }
      MPI_Bcast(&Rho_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /* READ IN OUTPUT FORMAT PARAMETERS */

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d  %d  %d %d ",&Lper_area, &Lcount_reflect, &Lprint_gofr, &Lprint_pmf);
    fprintf(fp2,"%d  %d  %d  %d",Lper_area,Lcount_reflect,Lprint_gofr,Lprint_pmf);
  }
  MPI_Bcast(&Lper_area,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Print_rho_type);
    fprintf(fp2,"%d",Print_rho_type);
  }
  MPI_Bcast(&Print_rho_type,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d %d",&Print_rho_switch,&Print_mesh_switch);
    fprintf(fp2,"%d  %d",Print_rho_switch,Print_mesh_switch);
  }
  MPI_Bcast(&Print_rho_type,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Iwrite);
    fprintf(fp2,"%d",Iwrite);
  }
  MPI_Bcast(&Iwrite,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) printf("\n TOTAL CHARGE IN dft_surfaces.dat = %9.6f\n",charge_sum);
  /* COARSENING Switches */

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Nzone);
    fprintf(fp2,"%d",Nzone);
  }
  MPI_Bcast(&Nzone,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nzone != 1) {
      for (izone =0; izone<Nzone-1; izone++){
        fscanf(fp,"%lf",&Rmax_zone[izone]);
        fprintf(fp2,"%f  ",Rmax_zone[izone]);
      if (Length_ref >0.0) Rmax_zone[izone] /= Length_ref;
      }
    }
  }
  MPI_Bcast(Rmax_zone,5,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Mesh_coarsening);
    fprintf(fp2,"%d",Mesh_coarsening);
  }
  MPI_Bcast(&Mesh_coarsening,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Mesh_coarsening==2 && Lauto_size && Proc==0){
      printf("AUTO SIZE FOR BOX :: [Lx,Ly,Lz]=");
      for (idim=0; idim<Ndim; idim++){
         Size_x[idim] = maxpos[idim]-minpos[idim] + 2.0*(Rmax_zone[0]+Sigma_ff[0][0]);
         printf("  %9.6f\n",Size_x[idim]);
      }
  }
  MPI_Bcast(Size_x,NDIM_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d %lf",&Coarser_jac,&Jac_grid);
    fprintf(fp2,"%d %f",Coarser_jac,Jac_grid);
  }
  MPI_Bcast(&Coarser_jac,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Jac_grid,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (Coarser_jac == 5) Nzone += 1;
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d %lf",&Lcut_jac,&Jac_threshold);
    fprintf(fp2,"%d %f",Lcut_jac,Jac_threshold);
  }
  MPI_Bcast(&Lcut_jac,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Jac_threshold,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /* Nonlinear Solver Parameters */

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d", &Max_Newton_iter);
    fprintf(fp2,"%d  ",Max_Newton_iter);
  }
  MPI_Bcast(&Max_Newton_iter,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%lg", &Newton_rel_tol);
    fprintf(fp2,"%lg  ",Newton_rel_tol);
    fscanf(fp,"%lg", &Newton_abs_tol);
    fprintf(fp2,"%lg  ",Newton_abs_tol);
    fscanf(fp,"%lg", &Min_update_frac);
    fprintf(fp2,"%lg  ",Min_update_frac);
  }
  MPI_Bcast(&Newton_rel_tol,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Newton_abs_tol,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Min_update_frac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d", &Load_Bal_Flag);
    fprintf(fp2,"%d  ", Load_Bal_Flag);
  }
  MPI_Bcast(&Load_Bal_Flag,1,MPI_INT,0,MPI_COMM_WORLD);

  /* Nonlinear Solver Parameters */

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d ", &L_Schur);
    fscanf(fp,"%d ", &Az_solver);
    if (Az_solver == 0) fscanf(fp,"%d ", &Az_kspace);
    else Az_kspace=-1;
    fprintf(fp2,"%d  %d ",Az_solver, Az_kspace);
  }
  MPI_Bcast(&L_Schur,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Az_solver,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Az_kspace,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d", &Az_scaling);
    fprintf(fp2,"%d  ",Az_scaling);
  }
  MPI_Bcast(&Az_scaling,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d  %lf", &Az_preconditioner, &Az_ilut_fill_param);
    fprintf(fp2,"%d  %f",Az_preconditioner, Az_ilut_fill_param);
  }
  MPI_Bcast(&Az_preconditioner,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Az_ilut_fill_param,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d  %lg", &Max_gmres_iter,&Az_tolerance);
    fprintf(fp2,"%d  %g  ",Max_gmres_iter,Az_tolerance);
  }
  MPI_Bcast(&Max_gmres_iter,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Az_tolerance,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /* READ IN MESH CONTINUATION PARAMETERS */

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Nruns);
    fprintf(fp2,"%d  ",Nruns);
  }
  MPI_Bcast(&Nruns,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Nruns > 0){

     if (Proc==0) {
       read_junk(fp,fp2);
       for (i=0; i < Ndim; ++i){
         fscanf(fp,"%lf", &Del_1[i]);
         fprintf(fp2,"%f  ",Del_1[i]);
       }
     }
     MPI_Bcast(Del_1,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

       if (Proc==0) {
         read_junk(fp,fp2);
         fscanf(fp,"%d",&Plane_new_nodes);
         fscanf(fp,"%d",&Pos_new_nodes);
         fprintf(fp2,"%d  %d",Plane_new_nodes,Pos_new_nodes);
       }
       MPI_Bcast(&Plane_new_nodes,1,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&Pos_new_nodes,1,MPI_INT,0,MPI_COMM_WORLD);
       if (Proc==0) {
         read_junk(fp,fp2);
           fscanf(fp,"%lf",&Guess_range[0]);
           fscanf(fp,"%lf",&Guess_range[1]);
	   fprintf(fp2,"%f  %f",Guess_range[0],Guess_range[1]);
       }
      MPI_Bcast(Guess_range,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
      
  }
  else{ read_junk(fp,fp2); read_junk(fp,fp2); read_junk(fp,fp2);}

    /* FINALLY READ IN LOCA PARAMETERS */
#ifdef USE_LOCA
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d", &itmp);
  }
  MPI_Bcast(&itmp,1,MPI_INT,0,MPI_COMM_WORLD);
  Loca.method = itmp;
  if (Proc==0) fprintf(fp2,"%d  ",Loca.method);
  if (Loca.method == 4) Lbinodal = 1;
  else Lbinodal=0;
  if (Lbinodal && Restart==0){
   printf("ERROR: Can only do binodal calculations with a restart\n") ;
   exit (-1);
  }

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d", &itmp);
    if (itmp == CONT_SCALE_RHO || itmp == CONT_SCALE_EPSW || itmp == CONT_SCALE_EPSWF ||
         itmp == CONT_SCALE_EPSFF || itmp== CONT_SCALE_CHG){
         fscanf(fp,"%lf ", &Scale_fac);
    }
    else Scale_fac = -1.;
  }
  MPI_Bcast(&itmp,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Scale_fac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  Loca.cont_type1 = itmp;
  if (Proc==0) fprintf(fp2,"%d  %g",Loca.cont_type1,Scale_fac);

  /* ALF: fix reading of input */
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%lf", &dtmp);
  }

  if (Proc==0 && Loca.method !=-1) {
    /* make adjustments for units */
    if (Loca.cont_type1 == CONT_TEMP && Temp>0.0 ) {
           dtmp /= (Eps_ff[0][0]*Temp);
    }
    else if ( Loca.cont_type1 == CONT_RHO_0 || Loca.cont_type1 == CONT_RHO_ALL ||
              Loca.cont_type1 == CONT_LOG_RHO_0 || Loca.cont_type1 == CONT_LOG_RHO_ALL ||  Loca.cont_type1 ==CONT_SCALE_RHO){
           if (Type_poly !=NONE && Type_poly !=WTC){
               printf("WARNING: density continuation may be incorrect with polymer code\n");
               printf("...... you really need a new c(r) for the hard chain part\n");
               printf("...... for each new density.\n");
           }
           if (Density_ref>0.0 && Loca.cont_type1 !=CONT_SCALE_RHO) dtmp /= Density_ref;
    }
    else if (Loca.cont_type1 == CONT_EPSW_0 ||  Loca.cont_type1 == CONT_EPSW_ALL ||
             Loca.cont_type1 == CONT_EPSWF00 ||  Loca.cont_type1 == CONT_EPSWF_ALL_0 ||
              Loca.cont_type1 == CONT_EPSFF_00 || Loca.cont_type1 == CONT_EPSFF_ALL ){
           if (Temp>0.0) dtmp /= Temp;
    }
    if (Mix_type == 0 && (Loca.cont_type1 == CONT_EPSWF00 ||  Loca.cont_type1 == CONT_EPSWF_ALL_0 ||
        Loca.cont_type1 == CONT_SCALE_EPSWF)){
        printf("ERROR: Can't do continuation in Eps_wf when the Mix_type is 0\n");
        exit(-1);
    } 
  }

    
  MPI_Bcast(&dtmp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  Loca.step_size = dtmp;
  if (Proc==0) fprintf(fp2,"%f  ",Loca.step_size);

  /* make adjustments to arrays if we are doing scaling continuation */
  if (Loca.cont_type1 ==  CONT_SCALE_RHO){
    for (icomp=0; icomp<Ncomp; icomp++) Rho_b[icomp] *= Scale_fac;
  }
  else if (Loca.cont_type1 ==  CONT_SCALE_EPSW){
    if (Mix_type==0) for (iwall_type=0; iwall_type<Nwall_type; iwall_type++) Eps_w[iwall_type] *= Scale_fac;
    else{
        for (i=0; i<Nwall_type;i++){
          for (j=0; j<Nwall_type; j++) Eps_ww[i][j] *= Scale_fac;
        }
    }
  }
  else if (Loca.cont_type1 ==  CONT_SCALE_EPSWF){
    for (icomp=0; icomp<Ncomp; icomp++){
       for (iwall_type=0; iwall_type<Nwall_type; iwall_type++) Eps_wf[icomp][iwall_type] *= Scale_fac;
    }
  }
  else if (Loca.cont_type1 ==  CONT_SCALE_EPSFF){
    if (Mix_type==0){ 
        for (icomp=0; icomp<Ncomp; icomp++) Eps_ff[icomp][icomp] *= Scale_fac;
    }
    else{
       for (i=0; i<Ncomp; i++){
         for (j=0; j<Ncomp;j++)  Eps_ff[i][j] *= Scale_fac;
       }
    }
  }
  else if (Loca.cont_type1 ==  CONT_SCALE_CHG){
    for (iwall=0; iwall<Nwall; iwall++) Elec_param_w[iwall] *= Scale_fac;
  }

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d", &itmp);
  }
  MPI_Bcast(&itmp,1,MPI_INT,0,MPI_COMM_WORLD);
  Loca.num_steps = itmp;
  if (Proc==0) fscanf(fp,"%lf", &dtmp);
  MPI_Bcast(&dtmp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  Loca.aggr = dtmp;
  if (Proc==0) fprintf(fp2,"%d %f ",Loca.num_steps, Loca.aggr);

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d", &itmp);
    if (itmp == CONT_SCALE_RHO || itmp == CONT_SCALE_EPSW || itmp == CONT_SCALE_EPSWF ||
         itmp == CONT_SCALE_EPSFF || itmp== CONT_SCALE_CHG){
         fscanf(fp,"%lf ", &Scale_fac);
    }
    else Scale_fac = -1.;
  }
  MPI_Bcast(&itmp,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Scale_fac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  Loca.cont_type2 = itmp;
  if (Proc==0) fprintf(fp2,"%d  %g",Loca.cont_type2,Scale_fac);

#else
  Loca.method = -1;
  if (Proc==0) fprintf(fp2,"Loca library not compiled in: %d  ",Loca.method);
#endif

  Nnodes_per_el_V = POW_INT(2, Ndim);
  Nnodes_per_el_S = POW_INT(2, Ndim-1);

  if (Proc==0) {
    fclose(fp);
    fclose(fp2);
  }
  error_check();

  /* with new units, we want to retain the ability to do continuation
     in temperature independent from the energy parameters Eps_ff or
     Eps_wf.  Therefore, we now set the reduced temperature for
     the purposes of continuation.  It is not used elsewhere in the
     calculations */
    
     if (Ipot_ff_n <= HARD_SPHERE) Temp = 1.0;
     else Temp = 1.0/Eps_ff[0][0];

  return;
}

/****************************************************************************/

/* read_junk :  This reads and then prints everything (i.e. the junk)
                between the last real input and the next.  It prints
                this junk to a file so that we will have saved the
                input file for the run of interest.                 */
void read_junk(FILE *fp, FILE *fp2)
{
   int c;
   while ((c=getc(fp)) != EOF && c !='@')
      if (Proc==0) putc(c,fp2);
}

/****************************************************************************/
/* error_check :  This function checks for error or inconsistencies in
                  the input parameters.  If either are found, the program
                  is stopped and error messages are sent to dft_errors.dat */

void error_check(void)
{
  double charge_b;
  int i,nmax;
  char *yo="error_check";

  if (Ndim> NDIM_MAX) {
    printf("%s: ERROR: Ndim (%d) larger than NDIM_MAX ", yo, Ndim);
    printf("(%d).\n\t\t Recompile with bigger NDIM_MAX\n", NDIM_MAX);
    exit(-1);
  }
  if (Nwall > NWALL_MAX){
     printf ("\n Nwall (%d) too large ... increase NWALL_MAX (%d)\n",Nwall,NWALL_MAX);
     exit(-1);
  }
  if (Nwall_type > NWALL_MAX_TYPE){
     printf ("\n Nwall_type (%d) too large ... increase NWALL_MAX_TYPE (%d)\n",Nwall,NWALL_MAX_TYPE);
     exit(-1);
  }
  if (Ncomp> NCOMP_MAX) {
    printf("%s: ERROR: Ncomp (%d) larger than NCOMP_MAX ", yo, Ncomp);
    printf("(%d).\n\t\t  Recompile with bigger NCOMP_MAX\n", NCOMP_MAX);
    exit(-1);
  }



  if (Ipot_ff_n > 2){
     printf ("\nSorry, your choice for the fluid-fluid interaction is not available\n");
     printf ("Ipot_ff_n: %d\n", Ipot_ff_n);
     exit (-1);
  }

  for (i=0;i<Nwall_type;i++){
     if (Ipot_wf_n[i] > 6 ){
        printf ("\nSorry, your choice for the wall-fluid interaction is not available\n");
        printf ("Ipot_wf_n[%d]: %d\n", i,Ipot_wf_n[i]);
        exit (-1);
     }
  }

  if ((Type_func > 2) || (Type_func < -1)){
     printf ("\nSorry, your choice for the hs functional is not yet available\n");
     printf ("Type_func: %d\n", Type_func);
     printf ("Reset Type_func to 0 for the Rosenfeld functional\n");
     exit (-1);
  }

  if (Type_attr > 0 && Type_attr < -1){
     printf ("\nSorry, your choice for the type of attractions is not yet available\n");
     printf ("Type_attr: %d\n", Type_attr);
     printf ("Reset Type_attr to 0 for the strict mean field approximation\n");
     exit (-1);
  }

  if (Type_coul > 1 && Type_coul < -1){
     printf ("\nSorry, your choice for the Coulomb functionals\n");
     printf ("Type_coul: %d\n", Type_coul);
     exit (-1);
  }

  if (Type_poly > WTC || Type_poly < NONE){
     printf ("\nSorry, your choice for the type of polymer functional\n");
     printf ("Type_poly: %d is not available\n", Type_poly);
     exit (-1);
  }

  charge_b = 0.0;
  if (Ipot_ff_c !=0) {
     for (i=0; i<Ncomp; i++) {
             charge_b += Rho_b[i]*Charge_f[i];
     }
     if (charge_b > 1.0e-10) {
        printf ("ERROR:input file: the fluid is not neutral\n");
        printf ("charge_b: %lf\n",charge_b);
        exit (-1);
     }
  }

  if (Ipot_ff_n == LJ12_6 && !Sten_Type[U_ATTRACT] && !Sten_Type[POLYMER_CR]) {
      printf("\nERROR: the u_attract func Sten_Type cannot be turned off\n");
      printf("if you wish to do an attracting fluid. Ipot_ff_n: %d\n", Ipot_ff_n);
      exit (-1);
  } 

  if (Ipot_ff_n != LJ12_6 && Sten_Type[U_ATTRACT] != FALSE) {
      printf("\nERROR: the u_attract func Sten_Type must be turned off\n");
      printf("if you don't want attractions in your calculation\n");
      printf("Ipot_ff_n=%d  Sten_Type[UATTRACT]=%d\n",Ipot_ff_n,Sten_Type[U_ATTRACT]);
      exit (-1);
  } 

  if (Load_Bal_Flag < 0 || Load_Bal_Flag > 5) {
      printf("\nERROR: the Load_Bal_Flag (%d) must 0-5 \n", Load_Bal_Flag);
      exit (-1);
  }
  if (Num_Proc == 1) Load_Bal_Flag = 0;

  if (Nzone > NZONE_MAX || Nzone <1 ){
     printf("\nERROR: Nzone out of range: Minimum val=1; Maximum val=%d\n",NZONE_MAX);
     exit(-1);
  }
  if (Coarser_jac == 5) nmax= Nzone-1;
  else                  nmax = Nzone;

  if ( Iliq_vap > -1 ) { 
    if (Ncomp > 1){
       printf("\nWarning: Can only do liquid-vapor equilibria for one species....will turn off coexistence calcs\n");
    }
    if (Ipot_ff_n == HARD_SPHERE){
       printf("\nERROR: There is no liquid-vapor transition for hard spheres\n");
       exit(-1);
    }
    if (Ipot_ff_n == IDEAL_GAS){
       printf("\nERROR: There is no liquid-vapor transition for the ideal gas\n");
       exit(-1);
    }
    if (Iliq_vap == 3 && Iguess1 != 6 ){
       printf("\nERROR: If Iliq_vap = 3 then Iguess should be 6\n");
       exit(-1);
    }
  }

  /* ALF: temporary change for polymers */
  if (Iguess1 == 6 && Iliq_vap !=3 && Lsteady_state == FALSE && !Sten_Type[POLYMER_CR]){
    printf("\nERROR: If Iguess = 6 then Iliq_vap should be 3\n");
    exit(-1);
  }

  if (Ncomp > 1) {
    if(Iguess1 == CONST_RHO_L || Iguess1 == CONST_RHO_V) {
       printf("\nERROR: Can't do liquid-vapor equilibria for more than one species\n");
       printf("\t Change Iguess type to %d  currently Iguess=%d\n",CONST_RHO,Iguess1);
       exit(-1);
    }
    else if(Iguess1 == EXP_RHO_L || Iguess1 == EXP_RHO_V) {
       printf("\nERROR: Can't do liquid-vapor equilibria for more than one species\n");
       printf("\t Change Iguess type to %d\n",EXP_RHO);
       exit(-1);
    }
    else if(Iguess1 == CHOP_RHO_L || Iguess1 == CHOP_RHO_V) {
       printf("\nERROR: Can't do liquid-vapor equilibria for more than one species\n");
       printf("\t Change Iguess type to %d\n",CHOP_RHO);
       exit(-1);
    }
  }
}
