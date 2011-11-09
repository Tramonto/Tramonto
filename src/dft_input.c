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
 *  FILE: dft_input.c
 *
 *  This file contains routines that read in the input file and
 *  initializes some flags based on the input.
 */

#include "dft_input.h"

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
   
   FILE *fp=NULL, *fp2=NULL, *fp3=NULL;

   char *yo = "read_input_file";
   char poly_file[20];
   int icomp, jcomp, iwall,iwall_type, idim, ipol,
       i, izone, j, jwall,new_wall,logical,ncharge, seg, block[NCOMP_MAX][NBLOCK_MAX],
       block_type[NBLOCK_MAX],pol_number, nlink_chk,irand,irand_range,itmp,
       dim_tmp,Lauto_center,Lauto_size,jmin=0,jmax=0,
       lzeros,latoms,ltrues,jwall_type,seg_tot;
   double rho_tmp[NCOMP_MAX],dtmp,charge_sum,minpos[3],maxpos[3];
   int iblock,jblock,read_periodic,read_wedge,read_rough,read_linear,itype_poly,repeat_type;


  
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

  /* Open the Files */
  if (Proc==0) {
    if( (fp  = fopen(input_file,"r")) == NULL) {
      printf("Can't open file %s\n", input_file);
      exit(-1);
    }
    if( (fp2 = fopen(output_file1,"w+")) == NULL) {
      printf("Can't open file %s\n", output_file1);
      exit(-1);
    }
    fprintf(fp2,"test out the printing %s\n",output_file1);
  }

  LDeBroglie=FALSE;
  LBulk=FALSE;
  Type_interface=UNIFORM_INTERFACE;

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
    if ( ((Type_bc[idim][0]==PERIODIC)&&(Type_bc[idim][1]!=PERIODIC))
         || ((Type_bc[idim][0]!=PERIODIC)&&(Type_bc[idim][1]==PERIODIC)) ) {
         printf("%s: ERROR: Both BC in dimension %d must be periodic if one is\n",
              yo, idim);
         exit(-1); 
    }
  }


  /* hard sphere functionals */
  if ( Proc==0 ) {
    read_junk(fp,fp2);
    fscanf(fp,"%d %d",&Type_func,&Type_hsdiam);
    fprintf(fp2,"%d %d",Type_func,Type_hsdiam);
  }
  MPI_Bcast(&Type_func,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Type_func >3 || Type_func<-1){
    if (Proc==0) printf("ERROR Type_func out of range - should be -1,0,1,2 or 3\n");
    exit(-1);
  }
  MPI_Bcast(&Type_hsdiam,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Type_hsdiam >2 || Type_hsdiam<0){
    if (Proc==0) printf("ERROR Type_hsdiam out of range - should be 0,1,2\n");
    exit(-1);
  }

  /* attractive functionals */
  if ( Proc==0 ) {
    read_junk(fp,fp2);
    fscanf(fp,"%d  %d",&Type_attr,&Type_pairPot);
    fprintf(fp2,"%d  %d",Type_attr,Type_pairPot);
  }
  MPI_Bcast(&Type_attr,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Type_pairPot,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Type_attr >6 || Type_attr<-1){
     if (Proc==0) printf("ERROR Type_attr=%d out of range - should be -1 <Type_attr<6\n",Type_attr);
     exit(-1);
  }

  /* coulombic functionals */
  if ( Proc==0 ) {
    read_junk(fp,fp2);
    fscanf(fp,"%d",&Type_coul);
    fprintf(fp2,"%d",Type_coul);
  }
  MPI_Bcast(&Type_coul,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Type_coul >4 || Type_coul<-1){
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

  if (Type_poly >WJDC3 || Type_poly<NONE){
     if (Proc==0) printf("ERROR Type_poly out of range (bounds are %d,%d)\n",NONE,WJDC3);
     exit(-1);
  }

  /* Read in or set if known the Potential Type Paramters */
  if (Type_func == -1 && (Type_poly==NONE || Type_poly==WTC)  ) Ipot_ff_n = IDEAL_GAS;
  else if (Type_attr == -1)                                     Ipot_ff_n = HARD_SPHERE;
  else if (Type_attr != NONE)      Ipot_ff_n = LJ12_6;
  else {
     printf("ERROR WITH Type_func and Type_attr selections and conversion to Ipot_ff_n parameter \n");
     exit (-1);
  }
  MPI_Bcast(&Ipot_ff_n,1,MPI_INT,0,MPI_COMM_WORLD);

 
  if (Type_coul >=0 || Type_pairPot==PAIR_COULOMB_CS || Type_pairPot==PAIR_COULOMB)  Ipot_ff_c =1; /* Coulombic Fluid */
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

 /*
  * Allocate the SurfaceGeometry variable to be a 3D array of structures
  * So it can be accessed as SGeom[iwall_type].param
  */

  SGeom = (struct SurfaceGeom_Struct *) array_alloc(1,Nwall_type, sizeof(struct SurfaceGeom_Struct));

  if (Nwall>0) Xtest_reflect_TF = (int **) array_alloc (2, Nlink,Ndim, sizeof(int));
  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall > 0){
      for (i=0; i < Nlink; i++)
	for (idim=0; idim< Ndim; idim++){
	  fscanf(fp,"%d",&Xtest_reflect_TF[i][idim]);
	  fprintf(fp2,"%d  ",Xtest_reflect_TF[i][idim]);
	}
    }
    else fprintf(fp2,"Xtest_reflect_TF n/a");
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
    else fprintf(fp2,"Surface_type n/a");
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
    else fprintf(fp2,"Orientation n/a");
  }
  if (Nwall_type > 0) 
    MPI_Bcast(Orientation,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);


  charge_sum=0.0;
  #ifndef _MSC_VER
    srandom(135649);
  #else
    srand(135649);
  #endif
  if (Nwall_type != 0){
    if ( Proc==0) {
      for (idim=0; idim<Ndim; idim++){ minpos[idim] = 1000.; maxpos[idim]=-1000.;}
      if( (fp3  = fopen("dft_surfaces.dat","r")) == NULL) {
	printf("Can't open file dft_surfaces.dat\n");
	exit(-1);
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
             #ifndef _MSC_VER
               irand = random();
             #else
               irand = rand();
             #endif
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
    else fprintf(fp2,"WallParam n/a");
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
    else fprintf(fp2,"WallParam_2 n/a");
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
    else fprintf(fp2,"WallParam_3 n/a");
  }
  if (Nwall_type > 0) 
    MPI_Bcast(WallParam_3,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (i=0; i < 3; ++i){
	fscanf(fp,"%d", &Lapply_offset[i]);
	fprintf(fp2,"%d  ",Lapply_offset[i]);
      }
    else fprintf(fp2,"Lapply_offset[] n/a");
  }
  if (Nwall_type > 0) 
    MPI_Bcast(Lapply_offset,3,MPI_INT,0,MPI_COMM_WORLD);

                 /* Surface Roughness */
  if (Proc==0) {
    read_junk(fp,fp2);
    read_rough=FALSE;
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%d", &Lrough_surf[iwall_type]);
	fprintf(fp2,"%d  ",Lrough_surf[iwall_type]);
        if (Lrough_surf[iwall_type]==TRUE) read_rough=TRUE;
      }
    else fprintf(fp2,"Lrough_surf n/a");
  }
  MPI_Bcast(&read_rough,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Nwall_type > 0) 
    MPI_Bcast(Lrough_surf,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0 && read_rough==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%lf", &Rough_param_max[iwall_type]);
	fprintf(fp2,"%f  ",Rough_param_max[iwall_type]);
      }
    else fprintf(fp2,"Rough_param_max n/a");
  }
  if (Nwall_type > 0 && read_rough==TRUE) 
    MPI_Bcast(Rough_param_max,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0 && read_rough==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%lf", &Rough_length[iwall_type]);
	fprintf(fp2,"%f  ",Rough_length[iwall_type]);
      }
    else fprintf(fp2,"Rough_length n/a");
  }
  if (Nwall_type > 0 && read_rough==TRUE) 
    MPI_Bcast(Rough_length,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

   /* populate a roughness array with random numbers that can be used to 
      generate roughness profiles */
  if (read_rough){
  for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){
      for (iblock=0;iblock<MAX_ROUGH_BLOCK;iblock++){
         for (jblock=0;jblock<MAX_ROUGH_BLOCK;jblock++){
            #ifndef _MSC_VER
              irand = random();
            #else
              irand = rand();
            #endif
            irand_range = POW_INT(2,31)-1;
            Rough_precalc[iwall_type][iblock][jblock]= Rough_param_max[iwall_type]*(-0.5+( ((double)irand)/((double)irand_range)));
         }
      }
  }
  }

                 /* Periodic Overlay  Params */
  if (Proc==0) {
    read_junk(fp,fp2);
    read_periodic=FALSE;
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
        Lperiodic_overlay[iwall_type]=FALSE;
	fscanf(fp,"%d", &Nperiodic_overlay[iwall_type]);
	fprintf(fp2,"%d  ",Nperiodic_overlay[iwall_type]);
        if (Nperiodic_overlay[iwall_type]>0){
            Lperiodic_overlay[iwall_type]=TRUE;
            read_periodic=TRUE;
        }
      }
    else fprintf(fp2,"Nperiodic_overlay n/a");
  }
  MPI_Bcast(&read_periodic,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Nwall_type > 0) {
    MPI_Bcast(Nperiodic_overlay,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Lperiodic_overlay,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
  }

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0&&read_periodic==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         for (i=0; i < Nperiodic_overlay[iwall_type]; ++i){
  	     fscanf(fp,"%d", &OrientationPeriodicFunc[iwall_type][i]);
     	     fprintf(fp2,"%d  ",OrientationPeriodicFunc[iwall_type][i]);
         }
      }
    else fprintf(fp2,"OrientationPeriodicFunc n/a");
  }
  if (Nwall_type > 0&&read_periodic==TRUE) 
    MPI_Bcast(OrientationPeriodicFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0&&read_periodic==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         for (i=0; i < Nperiodic_overlay[iwall_type]; ++i){
	    fscanf(fp,"%lf", &AmplitudePeriodicFunc[iwall_type][i]);
	    fprintf(fp2,"%f  ",AmplitudePeriodicFunc[iwall_type][i]);
         }
      }
    else fprintf(fp2,"AmplitudePeriodicFunc n/a");
  }
  if (Nwall_type > 0&&read_periodic==TRUE) 
    MPI_Bcast(AmplitudePeriodicFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0&&read_periodic==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         for (i=0; i < Nperiodic_overlay[iwall_type]; ++i){
	     fscanf(fp,"%lf", &WavelengthPeriodicFunc[iwall_type][i]);
	     fprintf(fp2,"%f  ",WavelengthPeriodicFunc[iwall_type][i]);
         }
      }
    else fprintf(fp2,"WavelengthPeriodicFunc n/a");
  }
  if (Nwall_type > 0&&read_periodic==TRUE) 
    MPI_Bcast(WavelengthPeriodicFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0&&read_periodic==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         for (i=0; i < Nperiodic_overlay[iwall_type]; ++i){
	     fscanf(fp,"%lf", &OriginPeriodicFunc[iwall_type][i]);
	     fprintf(fp2,"%f  ",OriginPeriodicFunc[iwall_type][i]);
         } 
      }
    else fprintf(fp2,"OrginPeriodicFunc n/a");
  }
  if (Nwall_type > 0&&read_periodic==TRUE) 
    MPI_Bcast(OriginPeriodicFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

                 /* Linear Overlay  Params */
  if (Proc==0) {
    read_junk(fp,fp2);
    read_linear=FALSE;
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
        Llinear_overlay[iwall_type]=FALSE;
	fscanf(fp,"%d", &Nlinear_overlay[iwall_type]);
	fprintf(fp2,"%d  ",Nlinear_overlay[iwall_type]);
        if (Nlinear_overlay[iwall_type]>0){
            Llinear_overlay[iwall_type]=TRUE;
            read_linear=TRUE;
        }
      }
    else fprintf(fp2,"Nlinear_overlay n/a");
  }
  MPI_Bcast(&read_linear,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Nwall_type > 0) {
    MPI_Bcast(Nlinear_overlay,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Llinear_overlay,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
  }

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0&&read_linear==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         for (i=0; i < Nlinear_overlay[iwall_type]; ++i){
  	     fscanf(fp,"%d", &OrientationLinearFunc[iwall_type][i]);
     	     fprintf(fp2,"%d  ",OrientationLinearFunc[iwall_type][i]);
         }
      }
    else fprintf(fp2,"OrientationLinearFunc n/a");
  }
  if (Nwall_type > 0&&read_linear==TRUE) 
    MPI_Bcast(OrientationLinearFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0&&read_linear==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         for (i=0; i < Nlinear_overlay[iwall_type]; ++i){
	    fscanf(fp,"%lf", &SlopeLinearFunc[iwall_type][i]);
	    fprintf(fp2,"%f  ",SlopeLinearFunc[iwall_type][i]);
         }
      }
    else fprintf(fp2,"SlopeLinerFunc n/a");
  }
  if (Nwall_type > 0&&read_linear==TRUE) 
    MPI_Bcast(SlopeLinearFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0&&read_linear==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         for (i=0; i < Nlinear_overlay[iwall_type]; ++i){
	     fscanf(fp,"%lf", &OriginLinearFunc[iwall_type][i]);
	     fprintf(fp2,"%f  ",OriginLinearFunc[iwall_type][i]);
         }
      }
    else fprintf(fp2,"OriginLinearFunc n/a");
  }
  if (Nwall_type > 0&&read_linear==TRUE) 
    MPI_Bcast(OriginLinearFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0&&read_linear==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         for (i=0; i < Nlinear_overlay[iwall_type]; ++i){
	     fscanf(fp,"%lf", &EndpointLinearFunc[iwall_type][i]);
	     fprintf(fp2,"%f  ",EndpointLinearFunc[iwall_type][i]);
         } 
      }
    else fprintf(fp2,"EndpointLinearFunc n/a");
  }
  if (Nwall_type > 0&&read_linear==TRUE) 
    MPI_Bcast(EndpointLinearFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

                 /* Angle Cutout Params */
  if (Proc==0) {
    read_junk(fp,fp2);
    read_wedge=FALSE;
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%d", &Lwedge_cutout[iwall_type]);
	fprintf(fp2,"%d  ",Lwedge_cutout[iwall_type]);
        if (Lwedge_cutout[iwall_type]==TRUE) read_wedge=TRUE;
      }
    else fprintf(fp2,"Lwedge_cutout n/a");
  }
  MPI_Bcast(&read_wedge,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Nwall_type > 0) 
    MPI_Bcast(Lwedge_cutout,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0&& read_wedge==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%lf", &Angle_wedge_start[iwall_type]);
	fprintf(fp2,"%f  ",Angle_wedge_start[iwall_type]);
      }
    else fprintf(fp2,"Angle_wedge_start n/a");
  }
  if (Nwall_type > 0&& read_wedge==TRUE) 
    MPI_Bcast(Angle_wedge_start,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0&& read_wedge==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fp,"%lf", &Angle_wedge_end[iwall_type]);
	fprintf(fp2,"%f  ",Angle_wedge_end[iwall_type]);
      }
    else fprintf(fp2,"Angle_wedge_end n/a");
  }
  if (Nwall_type > 0&& read_wedge==TRUE) 
    MPI_Bcast(Angle_wedge_end,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

             /* end surface adjustments */

  /* switches for types of wall-fluid and wall-wall interaction parameters */
  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         fscanf(fp,"%d",&Ipot_wf_n[iwall_type]);
         fprintf(fp2,"%d",Ipot_wf_n[iwall_type]);
      }
    else fprintf(fp2,"Ipot_wf_n n/a");
  }
  if (Nwall_type > 0) MPI_Bcast(Ipot_wf_n,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
 
          /* set logical indicating if any of the surfaces have hard cores - if so, we
              will need be careful with rosenfeld integrals */
  if (Proc==0){
     read_junk(fp,fp2);
     Lhard_surf=FALSE;
     if (Nwall_type>0){
        fscanf(fp,"%d",&Lhard_surf);
        fprintf(fp2,"%d ",Lhard_surf);
      }
      else fprintf(fp2,"Lhard_surf n/a");
  }
  if (Nwall_type>0){
      MPI_Bcast(&Lhard_surf,1,MPI_INT,0,MPI_COMM_WORLD);
  }

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         fscanf(fp,"%d",&Type_vext[iwall_type]);
         fprintf(fp2,"%d",Type_vext[iwall_type]);
      }
    else fprintf(fp2,"Type_vext n/a");
  }
  if (Nwall_type > 0) MPI_Bcast(Type_vext,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         fscanf(fp,"%d",&Vext_PotentialID[iwall_type]);
         fprintf(fp2,"%d",Vext_PotentialID[iwall_type]);
      }
    else fprintf(fp2,"Type_vext n/a");
  }
  if (Nwall_type > 0) MPI_Bcast(Vext_PotentialID,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    lzeros=FALSE; latoms=FALSE;
    if (Nwall_type > 0 && Ndim==3){ 
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
    }
    else fprintf(fp2,"Ipot_ww_n n/a");
  }
  MPI_Bcast(Ipot_ww_n,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0){
     read_junk(fp,fp2);
     if (Nwall_type > 0 && Ndim==3) {
        fscanf(fp,"%d",&Type_uwwPot);
        fprintf(fp2,"%d",Type_uwwPot);
     }
     else fprintf(fp2,"Type_uwwPot n/a");
  }
  if (Nwall_type > 0 && Ndim==3) MPI_Bcast(&Type_uwwPot,1,MPI_INT,0,MPI_COMM_WORLD);

  /* Fluid Particle Parameters */

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d  %d",&Ncomp,&Mix_type);
    fprintf(fp2,"%d  %d",Ncomp,Mix_type);
    if (Mix_type==0 && (Ipot_ff_n==HARD_SPHERE || Ipot_ff_n==IDEAL_GAS)){
        for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){
           if (Ipot_wf_n[iwall_type]!=VEXT_NONE && Ipot_wf_n[iwall_type]!=VEXT_HARD){
                printf("don't do LB mixing rules for hard sphere or ideal gas fluids if any\n");
                printf("of the surfaces has a nonzero Vext because Eps_f=0 and/or Sigma_f=0 by default\n");
                exit(-1);
           }
        }
    }
  }
  MPI_Bcast(&Ncomp,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Mix_type,1,MPI_INT,0,MPI_COMM_WORLD);

                      /* Manual HS_diam Entries */
  if (Proc==0) {
    read_junk(fp,fp2);
    if (Type_hsdiam==MANUAL_HS_DIAM){
       for (i=0; i<Ncomp; i++){
         fscanf(fp,"%lf",&HS_diam[i]); 
         fprintf(fp2,"%f  ",HS_diam[i]);
       }
    }
    else fprintf(fp2,"no read for HS_diam_manual_entry \n");
  }
  if (Type_hsdiam==MANUAL_HS_DIAM) MPI_Bcast(HS_diam,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

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
      fprintf(fp2,"Charge_f n/a");
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
      fprintf(fp2,"Sigma_ff n/a");
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

                       /* Yukawa parameters - the Debye length */
  if (Proc==0) {
    if (Type_pairPot == PAIR_YUKAWA_CS || Type_pairPot == PAIR_EXP_CS || 
        Type_pairPot==PAIR_LJandYUKAWA_CS || Type_pairPot==PAIR_r12andYUKAWA_CS ||
        Type_pairPot==PAIR_r18andYUKAWA_CS || Type_pairPot==PAIR_rNandYUKAWA_CS){
      read_junk(fp,fp2);
      for (i=0; i<Ncomp; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
        for (j=jmin; j<jmax; j++){
	  fscanf(fp,"%lf",&EpsYukawa_ff[i][j]);
	  fprintf(fp2,"%f  ",EpsYukawa_ff[i][j]);
          if (Temp > 0.0) EpsYukawa_ff[i][j]/=Temp;
        }
      }
      read_junk(fp,fp2);
      for (i=0; i<Ncomp; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
        for (j=jmin; j<jmax; j++){
	  fscanf(fp,"%lf",&YukawaK_ff[i][j]);
	  fprintf(fp2,"%f  ",YukawaK_ff[i][j]);
          if (Length_ref > 0.0) YukawaK_ff[i][j]*=Length_ref;
        }
      }
    }
    else  {
      read_junk(fp,fp2);
      for (i=0; i<Ncomp; i++){
         if (Mix_type==0) {jmin=i; jmax=i+1;}
         else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
         for(j=jmin;j<jmax;j++) EpsYukawa_ff[i][j] = 0.0;
      }
      fprintf(fp2,"EpsYukawa_ff n/a");
    
      read_junk(fp,fp2);
      for (i=0; i<Ncomp; i++){
         if (Mix_type==0) {jmin=i; jmax=i+1;}
         else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
         for(j=jmin;j<jmax;j++) YukawaK_ff[i][j] = 0.0;
      }
      fprintf(fp2,"YukawaK_ff n/a");
    }

    if(Type_pairPot==PAIR_rNandYUKAWA_CS){
      read_junk(fp,fp2);
      for (i=0; i<Ncomp; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
        for (j=jmin; j<jmax; j++){
	  fscanf(fp,"%lf",&Npow_ff[i][j]);
	  fprintf(fp2,"%f  ",Npow_ff[i][j]);
        }
      }
    }
    else{
      read_junk(fp,fp2);
      for (i=0; i<Ncomp; i++){
         if (Mix_type==0) {jmin=i; jmax=i+1;}
         else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
         for(j=jmin;j<jmax;j++) Npow_ff[i][j] = 0.0;
      }
      fprintf(fp2,"Npow_ff n/a");
    }
  }


  if (Type_pairPot==PAIR_YUKAWA_CS || Type_pairPot == PAIR_EXP_CS || 
      Type_pairPot==PAIR_LJandYUKAWA_CS || Type_pairPot==PAIR_r12andYUKAWA_CS ||
      Type_pairPot==PAIR_r18andYUKAWA_CS || Type_pairPot==PAIR_rNandYUKAWA_CS) {
          MPI_Bcast(EpsYukawa_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
	  MPI_Bcast(YukawaK_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
	  MPI_Bcast(Npow_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
      }

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

    /* resend WallParam to catch cases where it has been adjusted to be Sigma */
    MPI_Bcast(WallParam,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    fill_surfGeom_struct();

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

                       /* Yukawa parameters - the Debye length */
  if (Proc==0) {
    if (Type_uwwPot==PAIR_YUKAWA_CS || Type_uwwPot == PAIR_EXP_CS || 
        Type_uwwPot==PAIR_LJandYUKAWA_CS || Type_uwwPot==PAIR_r12andYUKAWA_CS
        || Type_uwwPot==PAIR_r18andYUKAWA_CS){

      read_junk(fp,fp2);
      for (i=0; i<Nwall_type; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Nwall_type;}
        for (j=jmin; j<jmax; j++){
	  if (Mix_type==1){ fscanf(fp,"%lf",&EpsYukawa_ww[i][j]);
	                    fprintf(fp2,"%f  ",EpsYukawa_ww[i][j]);
                            if (Temp > 0.0) EpsYukawa_ww[i][j]/=Temp;
                          }
          else            { fscanf(fp,"%lf",&EpsYukawa_w[i]);
	                    fprintf(fp2,"%f  ",EpsYukawa_w[i]);
                            if (Length_ref > 1.0) EpsYukawa_w[i]/=Temp;
                          }
        }
      }

      read_junk(fp,fp2);
      for (i=0; i<Nwall_type; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Nwall_type;}
        for (j=jmin; j<jmax; j++){
	  if (Mix_type==1){ fscanf(fp,"%lf",&YukawaK_ww[i][j]);
	                    fprintf(fp2,"%f  ",YukawaK_ww[i][j]);
                            if (Length_ref > 0.0) YukawaK_ww[i][j]*=Length_ref;
                          }
          else            { fscanf(fp,"%lf",&YukawaK_w[i]);
	                    fprintf(fp2,"%f  ",YukawaK_w[i]);
                            if (Length_ref > 0.0) YukawaK_w[i]*=Length_ref;
                          }
        }
      }
    }
    else{
      read_junk(fp,fp2);
      fprintf(fp2,"EpsYukawa_ww n/a");
      read_junk(fp,fp2);
      fprintf(fp2,"YukawaK_ww n/a");
    }
  }
  if (Type_uwwPot==PAIR_YUKAWA_CS || Type_uwwPot == PAIR_EXP_CS || 
      Type_uwwPot==PAIR_LJandYUKAWA_CS || Type_uwwPot==PAIR_r12andYUKAWA_CS
      || Type_uwwPot==PAIR_r18andYUKAWA_CS){
     if (Mix_type==1){
           MPI_Bcast(EpsYukawa_ww,NWALL_MAX_TYPE*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
           MPI_Bcast(YukawaK_ww,NWALL_MAX_TYPE*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
     }
     else{
          MPI_Bcast(EpsYukawa_w,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
          MPI_Bcast(YukawaK_w,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
     }
  }

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
                       /* Yukawa Energy paramter */
       read_junk(fp,fp2);
         for (i=0; i<Ncomp; i++){
           for (j=0; j<Nwall_type; j++){
   	  fscanf(fp,"%lf",&EpsYukawa_wf[i][j]);
   	  fprintf(fp2,"%f  ",EpsYukawa_wf[i][j]);
             if (Temp > 0.0) EpsYukawa_wf[i][j]/=Temp;
           }
         }
                       /* Yukawa paramter */
       read_junk(fp,fp2);
         for (i=0; i<Ncomp; i++){
           for (j=0; j<Nwall_type; j++){
   	  fscanf(fp,"%lf",&YukawaK_wf[i][j]);
   	  fprintf(fp2,"%f  ",YukawaK_wf[i][j]);
             if (Length_ref > 0.0) YukawaK_wf[i][j]*=Length_ref;
           }
         }
     }
    MPI_Bcast(Sigma_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Eps_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Cut_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(EpsYukawa_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(YukawaK_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  else{
    if (Proc==0){
         read_junk(fp,fp2);
         fprintf(fp2,"\n  USING L_B Mixing Rules --- WALL-FLUID INTERACTIONS COMPUTED BY CODE\n");
         fprintf(fp2,"........MANUAL INPUT DOES NOT APPLY \n");
         fprintf(fp2,"skipping this parameter  ");
         for (i=0; i<4; i++){
                read_junk(fp,fp2);
                fprintf(fp2,"skipping this parameter  ");
         }
    }
  }

  /* START POLYMER INPUT FOR ALL KINDS OF POLYMERS */

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
	  if (Proc==0) printf("Error: Must increase NBLOCK_MAX");
	  exit(-1);
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
	/*printf("Length (Nmer) of polymer %d = %d\n",pol_number,
	       Nmer[pol_number]);*/
      }
    }
    MPI_Bcast(Nmer,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
    if (Proc==0) {
      Ntype_mer = 0;
      read_junk(fp,fp2);
      seg_tot=0;
      for (i=0; i<NBLOCK_MAX; ++i){ Nmer_t_total[i]=0; Type_mer_to_Pol[i]=-1;}
      for (pol_number=0; pol_number<Npol_comp; ++pol_number){
        Poly_to_Ntype[pol_number]=0;
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
            if (Type_mer_to_Pol[block_type[i]]!=-1 && Type_mer_to_Pol[block_type[i]]!=pol_number){
               printf("Error: identical types on different chains are not allowed\n");
               exit(-1);
            }
            else Type_mer_to_Pol[block_type[i]]=pol_number;
            repeat_type=FALSE;
            for (itype_poly=0;itype_poly<Poly_to_Ntype[pol_number];itype_poly++) if (Poly_to_Type[pol_number][itype_poly]==block_type[i]) repeat_type=TRUE;
            if (repeat_type==FALSE){
                Poly_to_Type[pol_number][Poly_to_Ntype[pol_number]]=block_type[i];
                Poly_to_Ntype[pol_number]++;
            }
            seg++; seg_tot++;
          }
	  if (block_type[i] > Ntype_mer) Ntype_mer = block_type[i];
	}
      }
      Ntype_mer++;
      /*printf("Number of segment types (Ntype_mer) = %d\n",Ntype_mer);*/
    }
    MPI_Bcast(&Ntype_mer,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Nmer_t,NCOMP_MAX*NBLOCK_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Nmer_t_total,NBLOCK_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(SegChain2SegAll,NCOMP_MAX*NMER_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Type_mer,NCOMP_MAX*NMER_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Poly_to_Ntype,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Poly_to_Type,NCOMP_MAX*NBLOCK_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Type_mer_to_Pol,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
	  
	/* any grafted chains in system? */
	  /* first assume no */
	  for(pol_number=0; pol_number<Npol_comp; ++pol_number)
		  Grafted[pol_number] = 0;
	  if(Proc==0) {
		  read_junk(fp,fp2);
		  for (pol_number=0; pol_number<Npol_comp; ++pol_number){
			  fscanf(fp,"%d",&Grafted[pol_number]);
			  fprintf(fp2,"%d   ",Grafted[pol_number]);
		  }
		  read_junk(fp,fp2);
		  for (pol_number=0; pol_number<Npol_comp; ++pol_number){
			  if(Grafted[pol_number]) {
				fscanf(fp,"%d", &Graft_wall[pol_number]);
				fprintf(fp2,"%d  ",Graft_wall[pol_number]);
			  }
		  }
		  read_junk(fp,fp2);
		  for (pol_number=0; pol_number<Npol_comp; ++pol_number){
			  if(Grafted[pol_number]) {
				  fscanf(fp,"%lf",&Rho_g[pol_number]);
				  fprintf(fp2,"%f   ",Rho_g[pol_number]);
			  }
		  }
	  }
	  MPI_Bcast(Grafted,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
	   MPI_Bcast(Graft_wall,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
	  MPI_Bcast(Rho_g,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);

    if (Proc==0) {
      read_junk(fp,fp2);
      fscanf(fp,"%d", &Type_poly_arch);
      if (Type_poly_arch==POLY_ARCH_FILE && Type_poly != SCFT){
         fscanf(fp,"%s", poly_file);
         fprintf(fp2,"%s  ",poly_file);
      }
    }
    MPI_Bcast(&Type_poly_arch,1,MPI_INT,0,MPI_COMM_WORLD);

    /*now set up the chain architecture...*/
    setup_chain_architecture(poly_file,fp2);

      if (Type_poly == CMS){  /* this bit only applies to the CMS functional */
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
       }
       else{
          if (Proc==0) {
             read_junk(fp,fp2);
             fprintf(fp2,"\n NO LIQUID STATE INPUT FOR NON-CMS RUN\n");
             fprintf(fp2,"n/a   ");
             for (i=0; i<2; i++) {
                read_junk(fp,fp2);
                fprintf(fp2,"n/a   ");
             }       
          }       
       }
  }     /* END POLYMER INPUT FOR ALL KINDS OF POLYMERS */
  else{
    if (Proc==0) {
      read_junk(fp,fp2);
      fprintf(fp2,"\n POLYMER INPUT NOT RELEVENT FOR THIS RUN\n");
      fprintf(fp2,"not read   ");
      for (i=0; i<10; i++) {
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
        if (lzeros==FALSE  && ltrues==FALSE){
           fscanf(fp,"%d ", &Lsemiperm[i][j]);
           fprintf(fp2,"%d  ",Lsemiperm[i][j]);
        }
        if (i==0 && j==0){
           if (Lsemiperm[i][j]==-2) lzeros=TRUE;
           else if (Lsemiperm[i][j]==-1) ltrues=TRUE;
        }
        if (lzeros==TRUE) Lsemiperm[i][j]=0;
        else if (ltrues==TRUE) Lsemiperm[i][j]=1;
      }
    }
  }
  if (Nwall_type >0) MPI_Bcast(*Lsemiperm,Nwall_type*Ncomp,MPI_INT,0,MPI_COMM_WORLD);


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
    fscanf(fp,"%d  %d %d", &Type_interface, &Grad_dim,&Lconstrain_interface);
    fprintf(fp2,"%d  %d  %d",Type_interface,Grad_dim,Lconstrain_interface);
  }
  MPI_Bcast(&Type_interface,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Grad_dim,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Lconstrain_interface,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Proc==0) {
    read_junk(fp,fp2);
    if (Type_poly==NONE || Ntype_mer == 1){
      for (icomp=0; icomp<Ncomp; ++icomp){
        fscanf(fp,"%lf", &Rho_b[icomp]);
        fprintf(fp2,"%f  ",Rho_b[icomp]);
        if (Density_ref > 0.0) Rho_b[icomp] /= Density_ref;
        if (Type_interface != UNIFORM_INTERFACE) Rho_b_LBB[icomp]=Rho_b[icomp];
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
        if (Type_interface != UNIFORM_INTERFACE) Rho_b_LBB[i]=Rho_b[i];
        fprintf(fp2,"%f  ",Rho_b[i]);
      }
      for (icomp=Ntype_mer; icomp<Ncomp; ++icomp){
        fprintf(fp2,"%f  ",Rho_b[icomp]);
      }
    }
  }
  MPI_Bcast(Rho_b,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (Type_interface != UNIFORM_INTERFACE) MPI_Bcast(Rho_b_LBB,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /* calculate total sum of site densities */
  Rho_t = 0.0;
  for(icomp=0; icomp<Ncomp; icomp++) Rho_t += Rho_b[icomp];
  /*if (Type_interface==UNIFORM_INTERFACE) printf("Rho_t = %f\n", Rho_t);
  else printf("Rho_t(left) = %f\n", Rho_t);*/

  if (Type_interface != UNIFORM_INTERFACE){
     if (Proc==0) {
       read_junk(fp,fp2);
       if (Type_poly==NONE || Ntype_mer == 1){
         for (icomp=0; icomp<Ncomp; ++icomp){
           fscanf(fp,"%lf", &Rho_b[icomp]);
           fprintf(fp2,"%f  ",Rho_b[icomp]);
           if (Density_ref > 0.0) Rho_b[icomp] /= Density_ref;
           if (Type_interface != UNIFORM_INTERFACE) Rho_b_RTF[icomp]=Rho_b[icomp];
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
           if (Type_interface != UNIFORM_INTERFACE) Rho_b_RTF[i]=Rho_b[i];
           fprintf(fp2,"%f  ",Rho_b[i]);
         }
         for (icomp=Ntype_mer; icomp<Ncomp; ++icomp){
           fprintf(fp2,"%f  ",Rho_b[icomp]);
         }
       }
     }
     MPI_Bcast(Rho_b_RTF,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
     /* calculate total sum of site densities */
     Rho_t = 0.0;
     for(icomp=0; icomp<Ncomp; icomp++) Rho_t += Rho_b_RTF[icomp];
	  printf("Rho_t(right) = %f\n", Rho_t);
  }
  else{
    if (Proc==0) {
      read_junk(fp,fp2);
      fprintf(fp2,"n/a");
    }
  }

  if (Ipot_ff_c == COULOMB && Type_interface != UNIFORM_INTERFACE) {
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

  if (Type_interface != UNIFORM_INTERFACE){
     if (Proc==0) {
         read_junk(fp,fp2);
         fscanf(fp,"%lf", &X_const_mu);
         fprintf(fp2,"%f  ",X_const_mu);
         if (Length_ref > 0.0) X_const_mu /= Length_ref;
     }
     MPI_Bcast(&X_const_mu,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  else{
      if (Proc==0) {
	read_junk(fp,fp2);
	fprintf(fp2,"n/a");
      }
  }


  /* Read in Charged surface parameters */
  Ipot_wf_c = 0;
  if (Type_coul != NONE || Type_pairPot==PAIR_COULOMB_CS || Type_pairPot==PAIR_COULOMB) {
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
      fscanf(fp,"%d  %lf %lf %lf ",&Type_dielec, &Sigma_Angstroms_plasma, &Temp_K_plasma, &DielecConst_plasma);
      fprintf(fp2,"%d  %f %f %f",Type_dielec,Sigma_Angstroms_plasma,Temp_K_plasma,DielecConst_plasma);
    }
    MPI_Bcast(&Type_dielec,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&Sigma_Angstroms_plasma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&Temp_K_plasma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&DielecConst_plasma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    Temp_elec = 4.0*PI*KBOLTZ*Temp_K_plasma*DielecConst_plasma*EPSILON_0*Sigma_Angstroms_plasma*1.e-10/(E_CONST*E_CONST);

       /*Temp_elec = 4*PI*KBOLTZ*298.0*KAPPA_H2O*EPSILON_0*4.25e-10/(E_CONST*E_CONST);  Tang-Davis Paper Parameters*/

    if (Proc==0) printf("\t plasma parameter=%9.6f\n",1./Temp_elec);

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
  
  /* READ IN DIFFUSION SPECIFIC PARAMETERS */
  L1D_bc = FALSE;
  Linear_transport=0;   /* only want generalized Fick's Law */
                        /* can do J=-Dgrad mu, but doesn't work well */
  MPI_Bcast(&Linear_transport,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Type_interface==DIFFUSIVE_INTERFACE) {

    if (Proc==0) {
      read_junk(fp,fp2);
      for (icomp=0; icomp<Ncomp; ++icomp){
	fscanf(fp,"%lf", &D_coef[icomp]);
	fprintf(fp2,"%f  ",D_coef[icomp]);
      }
    }
    MPI_Bcast(D_coef,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if (Proc==0) {
      read_junk(fp,fp2);
      fscanf(fp,"%lf ", &Velocity);
      fprintf(fp2,"%f  ",Velocity);
    }
    MPI_Bcast(&Velocity,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

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

  }
  else{
    if (Proc==0) {
      read_junk(fp,fp2);
      fprintf(fp2,"TRANSPORT INPUT NOT RELEVENT FOR THIS RUN\n");
      fprintf(fp2,"n/a");
      for (i=0; i<3; i++) {
	read_junk(fp,fp2);
	fprintf(fp2,"n/a");
      }
    }
  }

  /* ***********end of physics input .... now run control and numerics input ********/ 

  /* Run control parameters: initial guess type etc.*/

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d %d",&Iguess,&Iguess_fields);
    fprintf(fp2,"%d  %d ",Iguess,Iguess_fields);
  }
  MPI_Bcast(&Iguess,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Iguess_fields,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Iguess==STEP_PROFILE || (Iguess>=CHOP_RHO && Iguess<= CHOP_RHO_STEP)){
    if (Proc==0) {
      read_junk(fp,fp2);
      if (Iguess != STEP_PROFILE) Nsteps=1;
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
    fscanf(fp,"%d ",&Restart);
    if (Restart == RESTART_FEWERCOMP) fscanf(fp," %d",&Nmissing_densities);
    else Nmissing_densities=0;
    fprintf(fp2,"%d  %d",Restart,Nmissing_densities);
  }
  MPI_Bcast(&Restart,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Nmissing_densities,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Proc==0){
    read_junk(fp,fp2);
    fscanf(fp,"%d ",&Restart_Vext);
    fprintf(fp2,"%d ",Restart_Vext);
    if (Restart_Vext != READ_VEXT_FALSE){
        fscanf(fp,"%s", Vext_file);
        fprintf(fp2,"  %s ", Vext_file);
        if (Restart_Vext == READ_VEXT_SUMTWO || Restart_Vext == READ_VEXT_STATIC){
           fscanf(fp,"%s", Vext_file2);
           fprintf(fp2,"  %s ", Vext_file2);
        }
    }
  }
  MPI_Bcast(&Restart_Vext,1,MPI_INT,0,MPI_COMM_WORLD);
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
    if (Print_rho_type==0 && Restart != 0){
      printf("WARNING: Print_rho_type is being set to 1 so that restart files will not be overwritten\n");
      Print_rho_type=1;
    }
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

  if (Proc==0 && fabs(charge_sum) > 1.e-8 && Iwrite != NO_SCREEN) printf("\n TOTAL CHARGE IN dft_surfaces.dat = %9.6f\n",charge_sum);
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
      if(Iwrite != NO_SCREEN) printf("AUTO SIZE FOR BOX :: [Lx,Ly,Lz]=");
      for (idim=0; idim<Ndim; idim++){
         Size_x[idim] = maxpos[idim]-minpos[idim] + 2.0*(Rmax_zone[0]+Sigma_ff[0][0]);
         if(Iwrite != NO_SCREEN) printf("  %9.6f\n",Size_x[idim]);
      }
  }
  MPI_Bcast(Size_x,NDIM_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d %lf",&Coarser_jac,&Jac_grid);
    fprintf(fp2,"%d %f ",Coarser_jac,Jac_grid);
  }
  MPI_Bcast(&Coarser_jac,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Jac_grid,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (Coarser_jac == 5) Nzone += 1;
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d %lf",&Lcut_jac,&Jac_threshold);
    fprintf(fp2,"%d %f ",Lcut_jac,Jac_threshold);
  }
  MPI_Bcast(&Lcut_jac,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Jac_threshold,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (Proc==0) {
     read_junk(fp,fp2);
     fscanf(fp,"%d  %lf", &L1D_bc,&X_1D_bc );
     fprintf(fp2,"%d  %f  ",L1D_bc,X_1D_bc);
     if (Length_ref > 0.0) X_1D_bc /= Length_ref;
  }
  MPI_Bcast(&L1D_bc,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&X_1D_bc,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /* Nonlinear Solver Parameters */

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d %d %d %d %d", &NL_Solver, &Max_NL_iter, &Physics_scaling, &ATTInA22Block, &Analyt_WJDC_Jac);
    fprintf(fp2,"NL_Solver=%d %d %d %d %d",NL_Solver,Max_NL_iter,Physics_scaling,ATTInA22Block,Analyt_WJDC_Jac);
    read_junk(fp,fp2);
    if (Physics_scaling != FALSE  && (Type_poly == WJDC || Type_poly==WJDC2 || Type_poly==WJDC3)){
      for (ipol=0;ipol<Npol_comp;ipol++){
            for (icomp=0;icomp<Ncomp;icomp++) {
                  fscanf(fp,"%lf ", &Scale_fac_WJDC[ipol][icomp]);
                  fprintf(fp2,"%f ", Scale_fac_WJDC[ipol][icomp]);
            }
      }
    }
    else { if(Iwrite != NO_SCREEN) printf("n/a - no manual entry of scaling parameters"); }
  }
  MPI_Bcast(Scale_fac_WJDC,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Max_NL_iter,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&NL_Solver,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Physics_scaling,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&ATTInA22Block,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Analyt_WJDC_Jac,1,MPI_INT,0,MPI_COMM_WORLD);
  if (NL_Solver==PICARD_BUILT_IN && Iguess_fields !=CALC_ALL_FIELDS){
     if(Iwrite != NO_SCREEN) printf("Picard solver indicated so Iguess_fields is reset to %d\n",CALC_ALL_FIELDS);
  }

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%lg", &NL_rel_tol);
    fprintf(fp2,"%lg  ",NL_rel_tol);
    fscanf(fp,"%lg", &NL_abs_tol);
    fprintf(fp2,"%lg  ",NL_abs_tol);
    fscanf(fp,"%lg", &NL_update_scalingParam);
    fprintf(fp2,"%lg  ",NL_update_scalingParam);
  }
  MPI_Bcast(&NL_rel_tol,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&NL_abs_tol,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&NL_update_scalingParam,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%lg", &NL_rel_tol_picard);
    fprintf(fp2,"%lg  ",NL_rel_tol_picard);
    fscanf(fp,"%lg", &NL_abs_tol_picard);
    fprintf(fp2,"%lg  ",NL_abs_tol_picard);
  }
  MPI_Bcast(&NL_rel_tol_picard,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&NL_abs_tol_picard,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d", &Load_Bal_Flag);
    fprintf(fp2,"%d  ", Load_Bal_Flag);
  }
  MPI_Bcast(&Load_Bal_Flag,1,MPI_INT,0,MPI_COMM_WORLD);

  /* Linear Solver Parameters */

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d ", &L_Schur);
    fscanf(fp,"%d ", &Az_solver);
    if (Az_solver == 0) fscanf(fp,"%d ", &Az_kspace);
    else Az_kspace=-1;
    fprintf(fp2,"L_Schur=%d  %d  %d ", L_Schur, Az_solver, Az_kspace);
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
    fprintf(fp2,"%d  %f ",Az_preconditioner, Az_ilut_fill_param);
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
    fprintf(fp2,"Nruns=%d  ",Nruns);
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
         fprintf(fp2,"%d  %d ",Plane_new_nodes,Pos_new_nodes);
       }
       MPI_Bcast(&Plane_new_nodes,1,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&Pos_new_nodes,1,MPI_INT,0,MPI_COMM_WORLD);
  }
  else{ read_junk(fp,fp2); read_junk(fp,fp2);}

    /* FINALLY READ IN LOCA PARAMETERS */
#ifdef USE_LOCA
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d", &itmp);
    fprintf(fp2,"%d ",itmp);
  }
  MPI_Bcast(&itmp,1,MPI_INT,0,MPI_COMM_WORLD);
  Loca.method = itmp;
/*  if (Proc==0) fprintf(fp2,"%d  ",Loca.method);*/
  if (Loca.method == 4) Lbinodal = 1;
  else Lbinodal=0;
  if (Lbinodal && Restart==0){
   printf("ERROR: Can only do binodal calculations with a restart\n") ;
   exit (-1);
  }

  for (i=0;i<NCONT_MAX;i++){
     for (j=0;j<2;j++) Cont_ID[i][j]=-1;
  } 
  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%d %d ", &itmp, &NID_Cont);
    for (i=0;i<NID_Cont;i++) fscanf(fp,"%d ", &Cont_ID[0][i]);
  }
  MPI_Bcast(&itmp,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&NID_Cont,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(Cont_ID,NCONT_MAX*2,MPI_INT,0,MPI_COMM_WORLD);
  Loca.cont_type1 = itmp;
  if (Proc==0){
       fprintf(fp2,"%d  %d  ",Loca.cont_type1,NID_Cont);
       for (i=0;i<NID_Cont;i++) fprintf(fp2,"%d  ",Cont_ID[0][i]);
  }

  if (Proc==0) {
    read_junk(fp,fp2);
    fscanf(fp,"%lf", &dtmp);
  }

  if (Proc==0 && Loca.method !=-1) {
    /* make adjustments for units */
    if (Loca.cont_type1 == CONT_TEMP && Temp>0.0 ) {
           dtmp /= (Eps_ff[0][0]*Temp);
    }
    else if ( Loca.cont_type1 == CONT_RHO_I || 
              Loca.cont_type1 == CONT_LOG_RHO_I){
           if (Type_poly == CMS){
               printf("WARNING: density continuation may be incorrect with polymer code\n");
               printf("...... you really need a new c(r) for the hard chain part\n");
               printf("...... for each new density.\n");
           }
           if (Density_ref>0.0) dtmp /= Density_ref;
    }
    else if (Loca.cont_type1 == CONT_EPSW_I ||  Loca.cont_type1 == CONT_EPSW_ALL ||
             Loca.cont_type1 == CONT_EPSWF_IJ ||  Loca.cont_type1 == CONT_EPSWF_ALL ||
              Loca.cont_type1 == CONT_EPSFF_IJ || Loca.cont_type1 == CONT_EPSFF_ALL ){
           if (Temp>0.0) dtmp /= Temp;
    }
    if (Mix_type == 0 && (Loca.cont_type1 == CONT_EPSWF_IJ ||  
                          Loca.cont_type1 == CONT_EPSWF_ALL)){
        printf("ERROR: Can't do continuation in Eps_wf when the Mix_type is 0\n");
        exit(-1);
    } 
  }

  MPI_Bcast(&dtmp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  Loca.step_size = dtmp;
  if (Proc==0) fprintf(fp2,"%f  ",Loca.step_size);

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
    fscanf(fp,"%d %d ", &itmp, &NID_Cont);
    for (i=0;i<NID_Cont;i++) fscanf(fp,"%d ", &Cont_ID[1][i]);
  }
  MPI_Bcast(&itmp,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&NID_Cont,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(Cont_ID,NCONT_MAX*2,MPI_INT,0,MPI_COMM_WORLD);
  Loca.cont_type2 = itmp;
  if (Proc==0){
       fprintf(fp2,"%d  %d  ",Loca.cont_type2,NID_Cont);
       for (i=0;i<NID_Cont;i++) fprintf(fp2,"%d  ",Cont_ID[1][i]);
  }

  if (!LBulk && (((Loca.method != -1 && Loca.cont_type1 == CONT_BETAMU_I) ||
       (Loca.method==4 && Loca.cont_type2 == CONT_BETAMU_I)) 
      || ((Loca.method != -1 && Loca.cont_type1 == CONT_BETAMU_I_NEW) ||
          (Loca.method==4 && Loca.cont_type2 == CONT_BETAMU_I_NEW))) 
     ){
       /*printf("for continuation in chemical potential LBulk must be TRUE=%d .... resetting LBulk\n",Loca.cont_type1,1);*/
       LBulk=TRUE;
    }
  MPI_Bcast(&LBulk,1,MPI_INT,0,MPI_COMM_WORLD);

  /* checks on LBulk */
  /* first check that bulk boundaries are NOT used if LBulk=TRUE and chemical potentials will be varied */ 
  if (LBulk && Loca.method != -1 &&
     (Loca.cont_type1 == CONT_BETAMU_I) ){
     for (idim=0;idim<Ndim;idim++){
         for (i=0;i<2;i++){
            if (Type_bc[idim][i]==IN_BULK){ 
               if (Proc==0){
                   printf("Bulk boundary detected while LBulk=TRUE and continuation in chemical potential requested.\n");
                   printf("This will not work because Rho_b is not updated during a continuation run. \n");
                   printf("The boundary condition will be reset to REFLECT (2). \n");
               }
               Type_bc[idim][i]=LAST_NODE;
               Type_bc[idim][i]=REFLECT;
            }
         }
     } 
  }
  

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
     else{ Temp = 1.0/Eps_ff[0][0];}

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
   if (Proc==0) fprintf(fp2,"   junk:");
   while ((c=getc(fp)) != EOF && c !='@')
      if (Proc==0) putc(c,fp2);
   if (Proc==0) fprintf(fp2,"   data:");
}
/****************************************************************************/
/* fill_SurfGeom_struct : this function pushes all the data from the generic input
formats to the surface geometry structures that will be used to access data */
void fill_surfGeom_struct()
{
  int iw,idim,i;
  struct SurfaceGeom_Struct *sgeom_iw;
  double r,rsq;

  if (Nwall_type>0) Poly_graft_dist = (double *) array_alloc (1, Nwall_type, sizeof(double));

  for (iw=0;iw<Nwall_type;iw++){
    Poly_graft_dist[iw]=0.;
    sgeom_iw = &(SGeom[iw]); 
    sgeom_iw->surfaceTypeID=Surface_type[iw];
    sgeom_iw->Lperiodic_overlay=FALSE;
    sgeom_iw->Llinear_overlay=FALSE;
    sgeom_iw->Lwedge_cutout=FALSE;
    switch(sgeom_iw->surfaceTypeID)
    {
       case smooth_planar_wall:
            sgeom_iw->halfwidth = (double *) array_alloc(1, Ndim, sizeof(double));
            sgeom_iw->orientation=Orientation[iw];
            sgeom_iw->halfwidth[Orientation[iw]]=WallParam[iw];
            sgeom_iw->Lrough_surface=Lrough_surf[iw];
            if (Lrough_surf[iw]==TRUE){
               sgeom_iw->roughness=Rough_param_max[iw];
               sgeom_iw->roughness_length=Rough_length[iw];
            }
            sgeom_iw->Lwedge_cutout=Lwedge_cutout[iw];
            if(Lwedge_cutout[iw]==TRUE){
               sgeom_iw->angle_wedge_start=Angle_wedge_start[iw];
               sgeom_iw->angle_wedge_end=Angle_wedge_end[iw];
            }
            sgeom_iw->Lperiodic_overlay=Lperiodic_overlay[iw];
            sgeom_iw->Nperiodic_overlay=Nperiodic_overlay[iw];
            for (i=0;i<Nperiodic_overlay[iw]; i++){ 
               sgeom_iw->orientation_periodic[i]=OrientationPeriodicFunc[iw][i];
               sgeom_iw->amplitude[i]=AmplitudePeriodicFunc[iw][i];
               sgeom_iw->wavelength[i]=WavelengthPeriodicFunc[iw][i];
               sgeom_iw->origin_PeriodicFunc[i]=OriginPeriodicFunc[iw][i];
            }
            sgeom_iw->Llinear_overlay=Llinear_overlay[iw];
            sgeom_iw->Nlinear_overlay=Nlinear_overlay[iw];
            for (i=0;i<Nlinear_overlay[iw]; i++){ 
               sgeom_iw->orientation_linear[i]=OrientationLinearFunc[iw][i];
               sgeom_iw->slope[i]=SlopeLinearFunc[iw][i];
               sgeom_iw->origin_LinearFunc[i]=OriginLinearFunc[iw][i];
               sgeom_iw->endpoint_LinearFunc[i]=EndpointLinearFunc[iw][i];
            }
            Poly_graft_dist[iw]=sgeom_iw->halfwidth[sgeom_iw->orientation];
            break;

       case finite_planar_wall:
            sgeom_iw->halfwidth = (double *) array_alloc(1, Ndim, sizeof(double));
            sgeom_iw->halfwidth[0]=WallParam[iw];
            if (Ndim>1) sgeom_iw->halfwidth[1]=WallParam_2[iw];
            if (Ndim>2) sgeom_iw->halfwidth[2]=WallParam_3[iw];
            sgeom_iw->Lrough_surface=Lrough_surf[iw];
            if (Lrough_surf[iw]==TRUE){
               sgeom_iw->roughness=Rough_param_max[iw];
               sgeom_iw->roughness_length=Rough_length[iw];
            }
            sgeom_iw->Llinear_overlay=Llinear_overlay[iw];
            sgeom_iw->Nlinear_overlay=Nlinear_overlay[iw];
            for (i=0;i<Nlinear_overlay[iw]; i++){ 
               sgeom_iw->orientation_linear[i]=OrientationLinearFunc[iw][i];
               sgeom_iw->slope[i]=SlopeLinearFunc[iw][i];
               sgeom_iw->origin_LinearFunc[i]=OriginLinearFunc[iw][i];
               sgeom_iw->endpoint_LinearFunc[i]=EndpointLinearFunc[iw][i];
            }
            break;

       case point_surface:
            rsq=0.0;
            for (idim=0;idim<Ndim;idim++) rsq+=Esize_x[idim]*Esize_x[idim];
            r=sqrt(rsq)/2.0;
            sgeom_iw->radius=r;
            sgeom_iw->halfwidth = (double *) array_alloc(1, Ndim, sizeof(double));
            for(idim=0;idim<Ndim;idim++) sgeom_iw->halfwidth[idim]=0.5*Esize_x[idim];
            break;

       case colloids_cyl_sphere:
            sgeom_iw->radius=WallParam[iw];
            sgeom_iw->Lwedge_cutout=Lwedge_cutout[iw];
            if(Lwedge_cutout[iw]==TRUE){
               sgeom_iw->angle_wedge_start=Angle_wedge_start[iw];
               sgeom_iw->angle_wedge_end=Angle_wedge_end[iw];
            }
            sgeom_iw->Lrough_surface=Lrough_surf[iw];
            if (Lrough_surf[iw]==TRUE){
               sgeom_iw->roughness=Rough_param_max[iw];
               sgeom_iw->roughness_length=Rough_length[iw];
            }
            Poly_graft_dist[iw]=sgeom_iw->radius;
            break;

       case finite_cyl_3D:
            sgeom_iw->orientation=Orientation[iw];
            sgeom_iw->radius=WallParam[iw];
            sgeom_iw->halflength=WallParam_2[iw];
            sgeom_iw->Lwedge_cutout=Lwedge_cutout[iw];
            if(Lwedge_cutout[iw]==TRUE){
               sgeom_iw->angle_wedge_start=Angle_wedge_start[iw];
               sgeom_iw->angle_wedge_end=Angle_wedge_end[iw];
            }
            sgeom_iw->Lrough_surface=Lrough_surf[iw];
            if (Lrough_surf[iw]==TRUE){
               sgeom_iw->roughness=Rough_param_max[iw];
               sgeom_iw->roughness_length=Rough_length[iw];
            }
            sgeom_iw->Lperiodic_overlay=Lperiodic_overlay[iw];
            sgeom_iw->Nperiodic_overlay=Nperiodic_overlay[iw];
            for (i=0;i<Nperiodic_overlay[iw]; i++){ 
               sgeom_iw->orientation_periodic[i]=OrientationPeriodicFunc[iw][i];
               sgeom_iw->amplitude[i]=AmplitudePeriodicFunc[iw][i];
               sgeom_iw->wavelength[i]=WavelengthPeriodicFunc[iw][i];
               sgeom_iw->origin_PeriodicFunc[i]=OriginPeriodicFunc[iw][i];
               if (OrientationPeriodicFunc[iw][i] != Orientation[iw]){
                   printf("Orientation of periodic function must be the same as the orientation of the surface\n");
                   printf("for the 3D cylindrical surface.  Adjustments can only be made along length of cylinder\n");
                   printf("Resetting periodic orientation to match the surface orientation\n");
                   sgeom_iw->orientation_periodic[i]=Orientation[iw];
               }
            }
            sgeom_iw->Llinear_overlay=Llinear_overlay[iw];
            sgeom_iw->Nlinear_overlay=Nlinear_overlay[iw];
            for (i=0;i<Nlinear_overlay[iw]; i++){ 
               sgeom_iw->orientation_linear[i]=OrientationLinearFunc[iw][i];
               sgeom_iw->slope[i]=SlopeLinearFunc[iw][i];
               sgeom_iw->origin_LinearFunc[i]=OriginLinearFunc[iw][i];
               sgeom_iw->endpoint_LinearFunc[i]=EndpointLinearFunc[iw][i];
               if (OrientationLinearFunc[iw][i] != Orientation[iw]){
                   printf("Orientation of linear function must be the same as the orientation of the surface\n");
                   printf("for the 3D cylinder.  Adjustments can only be made along length of the surface.\n");
                   printf("Resetting linear orientation to match the surface orientation\n");
                   sgeom_iw->orientation_linear[i]=Orientation[iw];
               }
            }
            Poly_graft_dist[iw]=sgeom_iw->radius;
            break;

       case atomic_centers:
            break;

       case cyl2D_sphere3D_pore:
            sgeom_iw->radius=WallParam[iw];
            Poly_graft_dist[iw]=sgeom_iw->radius;
            sgeom_iw->Lrough_surface=Lrough_surf[iw];
            if (Lrough_surf[iw]==TRUE){
               sgeom_iw->roughness=Rough_param_max[iw];
               sgeom_iw->roughness_length=Rough_length[iw];
            }
            sgeom_iw->Lwedge_cutout=Lwedge_cutout[iw];
            if(Lwedge_cutout[iw]==TRUE){
               sgeom_iw->angle_wedge_start=Angle_wedge_start[iw];
               sgeom_iw->angle_wedge_end=Angle_wedge_end[iw];
            }
            break;

       case cyl3D_slit2D_pore:
            sgeom_iw->orientation=Orientation[iw];
            sgeom_iw->radius=WallParam[iw];
            sgeom_iw->halflength=WallParam_2[iw];
            sgeom_iw->Lrough_surface=Lrough_surf[iw];
            if (Lrough_surf[iw]==TRUE){
               sgeom_iw->roughness=Rough_param_max[iw];
               sgeom_iw->roughness_length=Rough_length[iw];
            }
            sgeom_iw->Lperiodic_overlay=Lperiodic_overlay[iw];
            sgeom_iw->Nperiodic_overlay=Nperiodic_overlay[iw];
            for (i=0;i<Nperiodic_overlay[iw]; i++){ 
               sgeom_iw->orientation_periodic[i]=OrientationPeriodicFunc[iw][i];
               sgeom_iw->amplitude[i]=AmplitudePeriodicFunc[iw][i];
               sgeom_iw->wavelength[i]=WavelengthPeriodicFunc[iw][i];
               sgeom_iw->origin_PeriodicFunc[i]=OriginPeriodicFunc[iw][i];
               if (OrientationPeriodicFunc[iw][i] != Orientation[iw]){
                   printf("Orientation of periodic function must be the same as the orientation of the surface\n");
                   printf("for the 3D cylindrical or 2D slit pore.  adjustments can only be made along length of pore\n");
                   printf("Resetting periodic orientation to match the surface orientation\n");
                   sgeom_iw->orientation_periodic[i]=Orientation[iw];
               }
            }
            sgeom_iw->Lwedge_cutout=Lwedge_cutout[iw];
            if(Lwedge_cutout[iw]==TRUE){
               sgeom_iw->angle_wedge_start=Angle_wedge_start[iw];
               sgeom_iw->angle_wedge_end=Angle_wedge_end[iw];
            }
            sgeom_iw->Llinear_overlay=Llinear_overlay[iw];
            sgeom_iw->Nlinear_overlay=Nlinear_overlay[iw];
            for (i=0;i<Nlinear_overlay[iw]; i++){ 
               sgeom_iw->orientation_linear[i]=OrientationLinearFunc[iw][i];
               sgeom_iw->slope[i]=SlopeLinearFunc[iw][i];
               sgeom_iw->origin_LinearFunc[i]=OriginLinearFunc[iw][i];
               sgeom_iw->endpoint_LinearFunc[i]=EndpointLinearFunc[iw][i];
               if (OrientationLinearFunc[iw][i] != Orientation[iw]){
                   printf("Orientation of linear function must be the same as the orientation of the surface\n");
                   printf("for the 3D cylindrical or 2D slit pore.  adjustments can only be made along length of pore\n");
                   printf("Resetting linear orientation to match the surface orientation\n");
                   sgeom_iw->orientation_linear[i]=Orientation[iw];
               }
            }
            Poly_graft_dist[iw]=sgeom_iw->radius;
            break;

       default: /* No surface found */
            printf("error with surface type iwall_type=%d not identified\n",iw);
            exit(-1);
            break;
    }
  }
  return;
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
     if (Ipot_wf_n[i] > 7 ){
        printf ("\nSorry, your choice for the wall-fluid interaction is not available\n");
        printf ("Ipot_wf_n[%d]: %d\n", i,Ipot_wf_n[i]);
        exit (-1);
     }
  }

  if ((Type_func > 3) || (Type_func < -1)){
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

  if (Type_coul > 4 && Type_coul < -1){
     printf ("\nSorry, your choice for the Coulomb functionals\n");
     printf ("Type_coul: %d is not available \n", Type_coul);
     exit (-1);
  }

  if (Type_poly > WJDC3 || Type_poly < NONE){
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

  if (Load_Bal_Flag < 0 || Load_Bal_Flag > 5) {
      printf("\nERROR: the Load_Bal_Flag (%d) must 0-5 \n", Load_Bal_Flag);
      exit (-1);
  }
  if (Num_Proc == 1) Load_Bal_Flag = 0;

  if (Nzone > NZONE_MAX || Nzone <1 ){
     printf("\nERROR: Nzone out of range: Minimum val=1; Maximum val=%d: Current val=%d\n",NZONE_MAX,Nzone);
     exit(-1);
  }
  if (Coarser_jac == 5) nmax= Nzone-1;
  else                  nmax = Nzone;

}
