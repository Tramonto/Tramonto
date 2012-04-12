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

void read_input_file(FILE *fpinput, FILE *fpecho)

/* Read the input file for Classical Fluids Density Functional Theory
 *
 *    Authors:  Laura Frink, 9225
 *              Andrew Salinger, 9221
 */

{
   /* Local variable declarations */
   
   FILE *fpsurfaces=NULL;

   char *yo = "read_input_file";
   int icomp, jcomp, iwall,iwall_type, idim, ipol,
       i, izone, j, jwall,new_wall,logical,ncharge, seg, block[NCOMP_MAX][NBLOCK_MAX],
       block_type[NBLOCK_MAX],pol_number, nlink_chk,irand,irand_range,itmp,
       dim_tmp,Lauto_center,Lauto_size,jmin=0,jmax=0,
       lzeros,latoms,ltrues,jwall_type,seg_tot;
   double rho_tmp[NCOMP_MAX],dtmp,charge_sum,minpos[3],maxpos[3];
   int iblock,jblock,itype_poly,repeat_type,graft_logical;
   char unk_char[20];
   char *densityFile_tmp,*densityFile2_tmp;


  
  /********************** BEGIN EXECUTION ************************************/

    printf("Reading Input File\n");
  
  /* Read in the Mesh, Surface, Potential Type, Fluid Particle, 
     Surface Particle, State Point, Functional, and Run Control Parameters.

     The Input file (dft_input.dat) may be freely formatted with comments
     as long as the @ symbol is placed at the beginning of all data input
     lines.

     The input file is copied to the output file dft_out.lis */

  /********************************************/
  /* Set a Few Constants                      */
  /********************************************/
  LDeBroglie=FALSE;
  LBulk=FALSE;
  Type_interface=UNIFORM_INTERFACE;

  /********************************************/
  /* Look for request for GUI                 */
  /********************************************/
   Open_GUI=FALSE;
   fgets(unk_char,5,fpinput);
   if (strncmp(unk_char,"GUI",3)==0){
     Open_GUI=TRUE;
   }

  /**************************************************************************/
  /* Set a directory for output and default density files for Tramonto run */
  /**************************************************************************/
/*      read_junk(fpinput,fpecho);
      fscanf(fpinput,"%s", OutputFileDir_array);
      fprintf(fpecho,"%s  ",OutputFileDir_array);
      OutputFileDir=OutputFileDir_array; */

      OutputFileDir=".";   /* just set to cwd for now...*/

  sprintf(DensityFile_array, "./dft_dens.dat");
  sprintf(DensityFile2_array, "./dft_dens2.dat");

  DensityFile=DensityFile_array;
  DensityFile2=DensityFile2_array;


  /********************************************/
  /* Initialize and Read Dimension Parameters */
  /********************************************/
  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%lf  %lf  %lf  %lf %lf", &Length_ref, &Density_ref, &Temp, &Dielec_ref, &VEXT_MAX);
   fprintf(fpecho,"%f  %f  %f  %f   %f", Length_ref,Density_ref,Temp,Dielec_ref,VEXT_MAX);

  /***************************************/
  /* Initialize and Read Mesh Parameters */
  /***************************************/

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d",&Ndim);
   fprintf(fpecho,"%d",Ndim);

  read_junk(fpinput,fpecho);
  for (idim=0; idim < Ndim; ++idim){
      fscanf(fpinput,"%lf", &Size_x[idim]);
       fprintf(fpecho,"%f  ",Size_x[idim]);
  }
  for (idim=Ndim; idim<NDIM_MAX; idim++) Size_x[idim]=0.0;

  read_junk(fpinput,fpecho);
  for (idim=0; idim < Ndim; ++idim){
    fscanf(fpinput,"%lf", &Esize_x[idim]);
     fprintf(fpecho,"%f  ",Esize_x[idim]);
  }
  for (idim=Ndim; idim<NDIM_MAX; idim++) Esize_x[idim]=0.0;

  for (idim=0; idim < Ndim; ++idim){
    read_junk(fpinput,fpecho);
    for (i=0; i < 2; i++){
	fscanf(fpinput,"%d", &Type_bc[idim][i]);
	 fprintf(fpecho,"%d  ",Type_bc[idim][i]);
    }
  }
                                         /* scroll through unneeded entries */
  if (Ndim == 1) { read_junk(fpinput,fpecho); read_junk(fpinput,fpecho); }
  else if (Ndim ==2) { read_junk(fpinput,fpecho); }

  for (idim=0; idim < Ndim; ++idim) {    /* error checking on Type_bc */
    if ( ((Type_bc[idim][0]==PERIODIC)&&(Type_bc[idim][1]!=PERIODIC))
         || ((Type_bc[idim][0]!=PERIODIC)&&(Type_bc[idim][1]==PERIODIC)) ) {
         printf("%s: ERROR: Both BC in dimension %d must be periodic if one is\n",
              yo, idim);
         exit(-1); 
    }
  }

  /*********************************************/
  /* Initialize and Read Functional Selections */
  /*********************************************/

  /***************************/
  /* hard sphere functionals */
  /***************************/
  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d %d",&Type_func,&Type_hsdiam);
   fprintf(fpecho,"%d %d",Type_func,Type_hsdiam);

  if (Type_func >3 || Type_func<-1){ printf("ERROR Type_func out of range - should be -1,0,1,2 or 3\n"); exit(-1); }
  if (Type_hsdiam >2 || Type_hsdiam<0){ printf("ERROR Type_hsdiam out of range - should be 0,1,2\n"); exit(-1); }

     /***************************/
     /* attractive functionals */
     /***************************/
  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d  %d",&Type_attr,&Type_pairPot);
   fprintf(fpecho,"%d  %d",Type_attr,Type_pairPot);

  if (Type_attr >6 || Type_attr<-1){ printf("ERROR Type_attr=%d out of range - should be -1 <Type_attr<6\n",Type_attr); exit(-1); }

     /*************************/
     /* coulombic functionals */
     /*************************/
  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d",&Type_coul);
   fprintf(fpecho,"%d",Type_coul);

  if (Type_coul >4 || Type_coul<-1){ printf("ERROR Type_coul out of range - should be -1,0,1\n"); exit(-1); }

     /***********************/
     /* polymer functionals */
     /***********************/
  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d",&Type_poly);
   fprintf(fpecho,"%d",Type_poly);

  if (Type_poly >WJDC3 || Type_poly<NONE){ printf("ERROR Type_poly out of range (bounds are %d,%d)\n",NONE,WJDC3); exit(-1); }

     /***********************************************************/
     /* Potential Type Paramters Related to Functional settings */
     /***********************************************************/
  if (Type_func == -1 && (Type_poly==NONE || Type_poly==WTC)  ) Ipot_ff_n = IDEAL_GAS;
  else if (Type_attr == -1)                                     Ipot_ff_n = HARD_SPHERE;
  else if (Type_attr != NONE)      Ipot_ff_n = LJ12_6;
  else {
     printf("ERROR WITH Type_func and Type_attr selections and conversion to Ipot_ff_n parameter \n");
     exit (-1);
  }
 
  if (Type_coul >=0 || Type_pairPot==PAIR_COULOMB_CS || Type_pairPot==PAIR_COULOMB)  Ipot_ff_c =1; /* Coulombic Fluid */
  else Ipot_ff_c=0;  /* Neutral Fluid */

  /*********************************************/
  /* Initialize and Read Surface Parameters    */
  /*********************************************/
  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d  %d  %d  %d  %d",&Nwall_type,&Nwall,&Nlink,&Lauto_center,&Lauto_size);
   fprintf(fpecho,"%d  %d  %d  %d  %d",Nwall_type,Nwall,Nlink,Lauto_center,Lauto_size);


  if (Nwall>0) Xtest_reflect_TF = (int **) array_alloc (2, Nlink,Ndim, sizeof(int));
  read_junk(fpinput,fpecho);
  if (Nwall > 0){
    for (i=0; i < Nlink; i++)
      for (idim=0; idim< Ndim; idim++){
         fscanf(fpinput,"%d",&Xtest_reflect_TF[i][idim]);
	  fprintf(fpecho,"%d  ",Xtest_reflect_TF[i][idim]);
      }
  }
  else  fprintf(fpecho,"Xtest_reflect_TF n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0) 
    for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
        fscanf(fpinput,"%d", &Surface_type[iwall_type]);
	 fprintf(fpecho,"%d  ",Surface_type[iwall_type]);
    }
  else  fprintf(fpecho,"Surface_type n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0) 
    for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fpinput,"%d", &Orientation[iwall_type]);
	 fprintf(fpecho,"%d  ",Orientation[iwall_type]);
    }
  else  fprintf(fpecho,"Orientation n/a");

                   /****************************/
                   /* set up surface positions */
                   /****************************/
  charge_sum=0.0;
  #ifndef _MSC_VER
    srandom(135649);
  #else
    srand(135649);
  #endif

  if (Nwall_type > 0){
    for (idim=0; idim<Ndim; idim++){ minpos[idim] = 1000.; maxpos[idim]=-1000.;}

    if( (fpsurfaces  = fopen("dft_surfaces.dat","r")) == NULL) {
	printf("Can't open file dft_surfaces.dat\n"); exit(-1);
    }

    for (iwall=0; iwall<Nwall; iwall++){
       fscanf(fpsurfaces,"%d  %d",&WallType[iwall], &Link[iwall]);
       for (idim=0; idim<Ndim; idim++) {
          dim_tmp=idim;
                      /* temporary rotation of coordinates 
                                if (idim==0) dim_tmp=1;
                                else if (idim==1) dim_tmp=2;
                                else if (idim==2) dim_tmp=0;*/
                      /* end of temporary code */
          fscanf(fpsurfaces,"%lf",&WallPos[dim_tmp][iwall]);

                                         /* below is code for random placement of surface */
                                         /* coordinates - it is accomplished with a known */
                                         /* flag value for the WallPos variable.  This needs */
                                         /* to be changed to a proper selection int he input file*/
          if (fabs(WallPos[dim_tmp][iwall]+9999.0)<1.e-6) { 
            #ifndef _MSC_VER
              irand = random();
            #else
              irand = rand();
            #endif
            irand_range = POW_INT(2,31)-1;
            WallPos[dim_tmp][iwall] = Size_x[idim]*(-0.5+( ((double)irand)/((double)irand_range)));
            printf("\n  Wall %d dim %d gets WallPos:%g \n",iwall,idim,WallPos[idim][iwall]);
             fprintf(fpecho,"\n Wall %d dim %d gets WallPos:%g \n",iwall,idim,WallPos[idim][iwall]);
          }   /* end random wall placement */


           if (WallPos[dim_tmp][iwall] < minpos[dim_tmp]) minpos[dim_tmp]=WallPos[dim_tmp][iwall];
           if (WallPos[dim_tmp][iwall] > maxpos[dim_tmp]) maxpos[dim_tmp]=WallPos[dim_tmp][iwall];
       }

       fscanf(fpsurfaces,"%lf",&Elec_param_w[iwall]);
       charge_sum+=Elec_param_w[iwall];
     } 

     /*for (idim=0; idim<Ndim; idim++) printf("\n idim: %d min pos: %9.6f max pos %9.6f \n",idim,minpos[idim],maxpos[idim]);*/
     if (Lauto_center){
       for (iwall=0; iwall<Nwall; iwall++) for (idim=0; idim<Ndim; idim++) WallPos[idim][iwall] -= 0.5*(maxpos[idim] + minpos[idim]);
     }
     fclose(fpsurfaces);

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
  }
  else{ Nwall = 0; Nlink = 0; }
  
                            /**************************************/
                            /* Set Up Surface Geometry Parameters */
                            /**************************************/

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0) 
     for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fpinput,"%lf", &WallParam[iwall_type]);
	 fprintf(fpecho,"%f  ",WallParam[iwall_type]);
        if (Surface_type[iwall_type]==point_surface) WallParam[iwall_type]=Esize_x[0];
  }
  else  fprintf(fpecho,"WallParam n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0) 
    for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
       fscanf(fpinput,"%lf", &WallParam_2[iwall_type]);
        fprintf(fpecho,"%f  ",WallParam_2[iwall_type]);
    }
  else  fprintf(fpecho,"WallParam_2 n/a");
  
  read_junk(fpinput,fpecho);
  if (Nwall_type > 0) 
    for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
       fscanf(fpinput,"%lf", &WallParam_3[iwall_type]);
        fprintf(fpecho,"%f  ",WallParam_3[iwall_type]);
    }
  else  fprintf(fpecho,"WallParam_3 n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0) 
    for (i=0; i < 3; ++i){
       fscanf(fpinput,"%d", &Lapply_offset[i]);
        fprintf(fpecho,"%d  ",Lapply_offset[i]);
    }
  else  fprintf(fpecho,"Lapply_offset[] n/a");

                 /* Surface Roughness */
  read_junk(fpinput,fpecho);
  read_rough=FALSE;
  if (Nwall_type > 0) 
    for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fpinput,"%d", &Lrough_surf[iwall_type]);
	 fprintf(fpecho,"%d  ",Lrough_surf[iwall_type]);
        if (Lrough_surf[iwall_type]==TRUE) read_rough=TRUE;
    }
  else  fprintf(fpecho,"Lrough_surf n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0 && read_rough==TRUE) 
     for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fpinput,"%lf", &Rough_param_max[iwall_type]);
	 fprintf(fpecho,"%f  ",Rough_param_max[iwall_type]);
    }
  else  fprintf(fpecho,"Rough_param_max n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0 && read_rough==TRUE) 
     for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fpinput,"%lf", &Rough_length[iwall_type]);
	 fprintf(fpecho,"%f  ",Rough_length[iwall_type]);
     }
  else  fprintf(fpecho,"Rough_length n/a");


                 /* Periodic Overlay  Params */
  read_junk(fpinput,fpecho);
  read_periodic=FALSE;
  if (Nwall_type > 0) 
    for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
        Lperiodic_overlay[iwall_type]=FALSE;
	fscanf(fpinput,"%d", &Nperiodic_overlay[iwall_type]);
         fprintf(fpecho,"%d  ",Nperiodic_overlay[iwall_type]);
        if (Nperiodic_overlay[iwall_type]>0){
          Lperiodic_overlay[iwall_type]=TRUE;
          read_periodic=TRUE;
         }
    }
  else  fprintf(fpecho,"Nperiodic_overlay n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0&&read_periodic==TRUE) 
    for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
       for (i=0; i < Nperiodic_overlay[iwall_type]; ++i){
           fscanf(fpinput,"%d", &OrientationPeriodicFunc[iwall_type][i]);
     	    fprintf(fpecho,"%d  ",OrientationPeriodicFunc[iwall_type][i]);
       }
    }
  else  fprintf(fpecho,"OrientationPeriodicFunc n/a");


  read_junk(fpinput,fpecho);
  if (Nwall_type > 0&&read_periodic==TRUE) 
    for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
       for (i=0; i < Nperiodic_overlay[iwall_type]; ++i){
	    fscanf(fpinput,"%lf", &AmplitudePeriodicFunc[iwall_type][i]);
	     fprintf(fpecho,"%f  ",AmplitudePeriodicFunc[iwall_type][i]);
       }
    }
  else  fprintf(fpecho,"AmplitudePeriodicFunc n/a");


  read_junk(fpinput,fpecho);
  if (Nwall_type > 0&&read_periodic==TRUE) 
    for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
       for (i=0; i < Nperiodic_overlay[iwall_type]; ++i){
	     fscanf(fpinput,"%lf", &WavelengthPeriodicFunc[iwall_type][i]);
	      fprintf(fpecho,"%f  ",WavelengthPeriodicFunc[iwall_type][i]);
       }
    }
  else  fprintf(fpecho,"WavelengthPeriodicFunc n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0&&read_periodic==TRUE) 
    for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         for (i=0; i < Nperiodic_overlay[iwall_type]; ++i){
	     fscanf(fpinput,"%lf", &OriginPeriodicFunc[iwall_type][i]);
	      fprintf(fpecho,"%f  ",OriginPeriodicFunc[iwall_type][i]);
     } 
  }
  else  fprintf(fpecho,"OrginPeriodicFunc n/a");


                 /* Linear Overlay  Params */
  read_junk(fpinput,fpecho);
  read_linear=FALSE;
  if (Nwall_type > 0) 
     for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
        Llinear_overlay[iwall_type]=FALSE;
	fscanf(fpinput,"%d", &Nlinear_overlay[iwall_type]);
	 fprintf(fpecho,"%d  ",Nlinear_overlay[iwall_type]);
        if (Nlinear_overlay[iwall_type]>0){
            Llinear_overlay[iwall_type]=TRUE;
            read_linear=TRUE;
        }
      }
  else  fprintf(fpecho,"Nlinear_overlay n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0&&read_linear==TRUE) 
     for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         for (i=0; i < Nlinear_overlay[iwall_type]; ++i){
  	     fscanf(fpinput,"%d", &OrientationLinearFunc[iwall_type][i]);
     	      fprintf(fpecho,"%d  ",OrientationLinearFunc[iwall_type][i]);
         }
     }
  else  fprintf(fpecho,"OrientationLinearFunc n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0&&read_linear==TRUE) 
     for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         for (i=0; i < Nlinear_overlay[iwall_type]; ++i){
	    fscanf(fpinput,"%lf", &SlopeLinearFunc[iwall_type][i]);
	     fprintf(fpecho,"%f  ",SlopeLinearFunc[iwall_type][i]);
         }
     }
  else  fprintf(fpecho,"SlopeLinerFunc n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0&&read_linear==TRUE) 
     for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
        for (i=0; i < Nlinear_overlay[iwall_type]; ++i){
	     fscanf(fpinput,"%lf", &OriginLinearFunc[iwall_type][i]);
	      fprintf(fpecho,"%f  ",OriginLinearFunc[iwall_type][i]);
         }
     }
  else  fprintf(fpecho,"OriginLinearFunc n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0&&read_linear==TRUE) 
     for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
        for (i=0; i < Nlinear_overlay[iwall_type]; ++i){
	     fscanf(fpinput,"%lf", &EndpointLinearFunc[iwall_type][i]);
	      fprintf(fpecho,"%f  ",EndpointLinearFunc[iwall_type][i]);
        } 
     }
  else  fprintf(fpecho,"EndpointLinearFunc n/a");


                 /* Angle Cutout Params */
  read_junk(fpinput,fpecho);
  read_wedge=FALSE;
  if (Nwall_type > 0) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	 fscanf(fpinput,"%d", &Lwedge_cutout[iwall_type]);
	  fprintf(fpecho,"%d  ",Lwedge_cutout[iwall_type]);
         if (Lwedge_cutout[iwall_type]==TRUE) read_wedge=TRUE;
       }
  else  fprintf(fpecho,"Lwedge_cutout n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0&& read_wedge==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fpinput,"%lf", &Angle_wedge_start[iwall_type]);
	 fprintf(fpecho,"%f  ",Angle_wedge_start[iwall_type]);
      }
  else  fprintf(fpecho,"Angle_wedge_start n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0&& read_wedge==TRUE) 
      for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
	fscanf(fpinput,"%lf", &Angle_wedge_end[iwall_type]);
	 fprintf(fpecho,"%f  ",Angle_wedge_end[iwall_type]);
      }
  else  fprintf(fpecho,"Angle_wedge_end n/a");


  /*******************************************************************************/
  /* Initialize and Read Wall-Fluid and Wall-Wall Interaction Type Parameters    */
  /*******************************************************************************/

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0) 
    for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         fscanf(fpinput,"%d",&Ipot_wf_n[iwall_type]);
          fprintf(fpecho,"%d",Ipot_wf_n[iwall_type]);
    }
  else  fprintf(fpecho,"Ipot_wf_n n/a");

          /* set logical indicating if any of the surfaces have hard cores - if so, we
              will need be careful with rosenfeld integrals */
  read_junk(fpinput,fpecho);
  Lhard_surf=FALSE;
  if (Nwall_type>0){
      fscanf(fpinput,"%d",&Lhard_surf);
       fprintf(fpecho,"%d ",Lhard_surf);
  }
  else  fprintf(fpecho,"Lhard_surf n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0) 
    for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         fscanf(fpinput,"%d",&Type_vext[iwall_type]);
          fprintf(fpecho,"%d",Type_vext[iwall_type]);
    }
  else  fprintf(fpecho,"Type_vext n/a");


  read_junk(fpinput,fpecho);
  if (Nwall_type > 0) 
     for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
         fscanf(fpinput,"%d",&Vext_PotentialID[iwall_type]);
          fprintf(fpecho,"%d",Vext_PotentialID[iwall_type]);
     }
  else  fprintf(fpecho,"Type_vext n/a");

  read_junk(fpinput,fpecho);
  lzeros=FALSE; latoms=FALSE;
  if (Nwall_type > 0 && Ndim==3){ 
     for (iwall_type=0; iwall_type < Nwall_type; ++iwall_type){
        for (jwall_type=0; jwall_type < Nwall_type;++jwall_type){
           if (!lzeros && !latoms){
             fscanf(fpinput,"%d",&Ipot_ww_n[iwall_type][jwall_type]);
              fprintf(fpecho,"%d",Ipot_ww_n[iwall_type][jwall_type]);
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
  else  fprintf(fpecho,"Ipot_ww_n n/a");

  read_junk(fpinput,fpecho);
  if (Nwall_type > 0 && Ndim==3) {
      fscanf(fpinput,"%d",&Type_uwwPot);
       fprintf(fpecho,"%d",Type_uwwPot);
  }
  else  fprintf(fpecho,"Type_uwwPot n/a");


  /**********************************************************/
  /* Initialize and Read Fluid-Fluid Interaction Parameters */
  /**********************************************************/

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d  %d",&Ncomp,&Mix_type);
   fprintf(fpecho,"%d  %d",Ncomp,Mix_type);
  if (Mix_type==0 && (Ipot_ff_n==HARD_SPHERE || Ipot_ff_n==IDEAL_GAS)){
      for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){
         if (Ipot_wf_n[iwall_type]!=VEXT_NONE && Ipot_wf_n[iwall_type]!=VEXT_HARD){
                printf("don't do LB mixing rules for hard sphere or ideal gas fluids if any\n");
                printf("of the surfaces has a nonzero Vext because Eps_f=0 and/or Sigma_f=0 by default\n");
                exit(-1);
         }
      }
  }

                      /* Manual HS_diam Entries */
  read_junk(fpinput,fpecho);
  if (Type_hsdiam==MANUAL_HS_DIAM){
     for (i=0; i<Ncomp; i++){
       fscanf(fpinput,"%lf",&HS_diam[i]); 
        fprintf(fpecho,"%f  ",HS_diam[i]);
     }
  }
  else  fprintf(fpecho,"no read for HS_diam_manual_entry \n");


                      /* MASS Entries */
  read_junk(fpinput,fpecho);
  for (i=0; i<Ncomp; i++){
    fscanf(fpinput,"%lf",&Mass[i]); 
     fprintf(fpecho,"%f  ",Mass[i]);
  }

                      /* VALENCE Entries */
  read_junk(fpinput,fpecho);
  if (Ipot_ff_c >0){
    for (icomp=0; icomp < Ncomp; ++icomp){
      fscanf(fpinput,"%lf", &Charge_f[icomp]);
       fprintf(fpecho,"%f  ",Charge_f[icomp]);
    }
  }
  else{ 
     fprintf(fpecho,"Charge_f n/a");
    for (icomp=0; icomp < Ncomp; ++icomp)Charge_f[icomp]=0.0;
  }

                      /* POLARIZATION Entries */
  read_junk(fpinput,fpecho);
  if (Type_coul ==2){
     for (icomp=0; icomp < Ncomp; ++icomp){
        fscanf(fpinput,"%lf", &Pol[icomp]);
         fprintf(fpecho,"%f  ",Pol[icomp]);
        if (Pol[icomp]!=0.0) Lpolarize[icomp]=TRUE;
        else                 Lpolarize[icomp]=FALSE;
     }
  }
  else  fprintf(fpecho,"n/a: no polarization");

                       /* Sigma_ff */
  read_junk(fpinput,fpecho);
  if (Ipot_ff_n != IDEAL_GAS){
    for (i=0; i<Ncomp; i++){
      if (Mix_type==0) {jmin=i; jmax=i+1;}
      else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
      for (j=jmin; j<jmax; j++){
           fscanf(fpinput,"%lf",&Sigma_ff[i][j]);
            fprintf(fpecho,"%f  ",Sigma_ff[i][j]);
      }
    }
  }
  else  {
    for (i=0; i<Ncomp; i++){
       if (Mix_type==0) {jmin=i; jmax=i+1;}
       else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
       for(j=jmin;j<jmax;j++) Sigma_ff[i][j] = 0.0;
    }
     fprintf(fpecho,"Sigma_ff n/a");
  }

                       /* Eps_ff */
  read_junk(fpinput,fpecho);
  if (Ipot_ff_n != IDEAL_GAS && Ipot_ff_n != HARD_SPHERE) {
    for (i=0; i<Ncomp; i++){
      if (Mix_type==0) {jmin=i; jmax=i+1;}
      else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
      for (j=jmin; j<jmax; j++){
            fscanf(fpinput,"%lf",&Eps_ff[i][j]);
             fprintf(fpecho,"%f  ",Eps_ff[i][j]);
      }
    }
  }
  else  {
    for (i=0; i<Ncomp; i++){ 
       if (Mix_type==0) {jmin=i; jmax=i+1;}
       else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
       for (j=jmin; j<jmax; j++) Eps_ff[i][j] = 0.0;
    }
     fprintf(fpecho,"n/a: no Eps_ff for ideal gas or hard sphere");
  }

                       /* Cut_ff */
  read_junk(fpinput,fpecho);
  if (Ipot_ff_n != IDEAL_GAS && Ipot_ff_n != HARD_SPHERE) {
    for (i=0; i<Ncomp; i++){
      if (Mix_type==0) {jmin=i; jmax=i+1;}
      else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
      for (j=jmin; j<jmax; j++){
         fscanf(fpinput,"%lf",&Cut_ff[i][j]);
          fprintf(fpecho,"%f  ",Cut_ff[i][j]);
      }
    }
  }
  else  {
    for (i=0; i<Ncomp; i++){ 
       if (Mix_type==0) {jmin=i; jmax=i+1;}
       else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
       for (j=jmin; j<jmax; j++) Cut_ff[i][j] = 0.0;
    }
     fprintf(fpecho,"n/a: no cutoff for ideal gas or hard sphere");
  }

                       /* Bond_ff */
  read_junk(fpinput,fpecho);
  if (Type_poly!=NONE){
    for (i=0; i<Ncomp; i++){
      if (Mix_type==0) {jmin=i; jmax=i+1;}
      else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
      for (j=jmin; j<jmax; j++){
         fscanf(fpinput,"%lf",&Bond_ff[i][j]);
          fprintf(fpecho,"%f  ",Bond_ff[i][j]);
      }
    }
  }
  else  {
    for (i=0; i<Ncomp; i++){
       if (Mix_type==0) {jmin=i; jmax=i+1;}
       else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
       for(j=jmin;j<jmax;j++) Bond_ff[i][j] = 0.0;
    }
     fprintf(fpecho,"n/a: bonds only needed for polymer runs");
  }

                       /* Yukawa parameters - the Debye length */
  if (Type_pairPot == PAIR_YUKAWA_CS || Type_pairPot == PAIR_EXP_CS || 
      Type_pairPot==PAIR_LJandYUKAWA_CS || Type_pairPot==PAIR_r12andYUKAWA_CS ||
      Type_pairPot==PAIR_r18andYUKAWA_CS || Type_pairPot==PAIR_rNandYUKAWA_CS){
      read_junk(fpinput,fpecho);
      for (i=0; i<Ncomp; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
        for (j=jmin; j<jmax; j++){
	  fscanf(fpinput,"%lf",&EpsYukawa_ff[i][j]);
	   fprintf(fpecho,"%f  ",EpsYukawa_ff[i][j]);
        }
      }
      read_junk(fpinput,fpecho);
      for (i=0; i<Ncomp; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
        for (j=jmin; j<jmax; j++){
	  fscanf(fpinput,"%lf",&YukawaK_ff[i][j]);
	   fprintf(fpecho,"%f  ",YukawaK_ff[i][j]);
        }
      }
  }
  else  {
      read_junk(fpinput,fpecho);
      for (i=0; i<Ncomp; i++){
         if (Mix_type==0) {jmin=i; jmax=i+1;}
         else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
         for(j=jmin;j<jmax;j++) EpsYukawa_ff[i][j] = 0.0;
      }
       fprintf(fpecho,"EpsYukawa_ff n/a");
    
      read_junk(fpinput,fpecho);
      for (i=0; i<Ncomp; i++){
         if (Mix_type==0) {jmin=i; jmax=i+1;}
         else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
         for(j=jmin;j<jmax;j++) YukawaK_ff[i][j] = 0.0;
      }
       fprintf(fpecho,"YukawaK_ff n/a");
  }

  if(Type_pairPot==PAIR_rNandYUKAWA_CS){
      read_junk(fpinput,fpecho);
      for (i=0; i<Ncomp; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
        for (j=jmin; j<jmax; j++){
	  fscanf(fpinput,"%lf",&Npow_ff[i][j]);
	   fprintf(fpecho,"%f  ",Npow_ff[i][j]);
        }
      }
  }
  else{
      read_junk(fpinput,fpecho);
      for (i=0; i<Ncomp; i++){
         if (Mix_type==0) {jmin=i; jmax=i+1;}
         else if (Mix_type==1) {jmin=0;jmax=Ncomp;}
         for(j=jmin;j<jmax;j++) Npow_ff[i][j] = 0.0;
      }
       fprintf(fpecho,"Npow_ff n/a");
  }

  /********************************************************/
  /* Initialize and Read Wall-Wall Interaction Parameters */
  /********************************************************/

                       /* Density of atoms in the surfaces */
  read_junk(fpinput,fpecho);
  for (i=0; i<Nwall_type; i++){
      fscanf(fpinput,"%lf", &Rho_w[i]);
       fprintf(fpecho,"%f  ",Rho_w[i]);
  }
    
                       /*Sigma_w or Sigma_ww*/
  read_junk(fpinput,fpecho);
  for (i=0; i<Nwall_type; i++){
      if (Mix_type==0) {jmin=i; jmax=i+1;}
      else if (Mix_type==1) {jmin=0;jmax=Nwall_type;}

      for (j=jmin; j<jmax; j++){
          if (Mix_type==1){ 
             fscanf(fpinput,"%lf",&Sigma_ww[i][j]);
              fprintf(fpecho,"%f  ",Sigma_ww[i][j]);
             if (Surface_type[i] == atomic_centers && j==i) WallParam[i] = Sigma_ww[i][i]/2.0;
          }
          else { 
             fscanf(fpinput,"%lf",&Sigma_w[i]);
              fprintf(fpecho,"%f  ",Sigma_w[i]);
             if (Surface_type[i] == atomic_centers) WallParam[i] = Sigma_w[i]/2.0;
          }
      }
  }

                       /*Eps_w or Eps_ww*/
  read_junk(fpinput,fpecho);
  for (i=0; i<Nwall_type; i++){
      if (Mix_type==0) {jmin=i; jmax=i+1;}
      else if (Mix_type==1) {jmin=0;jmax=Nwall_type;}

      for (j=jmin; j<jmax; j++){
          if (Mix_type==1){ 
              fscanf(fpinput,"%lf",&Eps_ww[i][j]);
	       fprintf(fpecho,"%f  ",Eps_ww[i][j]);
          }
          else { 
              fscanf(fpinput,"%lf",&Eps_w[i]);
	       fprintf(fpecho,"%f  ",Eps_w[i]);
          }
        }
  }

                       /*Cut_ww */
  read_junk(fpinput,fpecho);
  for (i=0; i<Nwall_type; i++){
      if (Mix_type==0) {jmin=i; jmax=i+1;}
      else if (Mix_type==1) {jmin=0;jmax=Nwall_type;}

      for (j=jmin; j<jmax; j++){
          fscanf(fpinput,"%lf",&Cut_ww[i][j]);
	   fprintf(fpecho,"%f  ",Cut_ww[i][j]);
      }
  }

                       /* Yukawa parameters - the Debye length */
  if (Type_uwwPot==PAIR_YUKAWA_CS || Type_uwwPot == PAIR_EXP_CS || 
      Type_uwwPot==PAIR_LJandYUKAWA_CS || Type_uwwPot==PAIR_r12andYUKAWA_CS
      || Type_uwwPot==PAIR_r18andYUKAWA_CS){

      read_junk(fpinput,fpecho);
      for (i=0; i<Nwall_type; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Nwall_type;}
        for (j=jmin; j<jmax; j++){
	  if (Mix_type==1){ fscanf(fpinput,"%lf",&EpsYukawa_ww[i][j]);
	                     fprintf(fpecho,"%f  ",EpsYukawa_ww[i][j]);
                          }
          else            { fscanf(fpinput,"%lf",&EpsYukawa_w[i]);
	                     fprintf(fpecho,"%f  ",EpsYukawa_w[i]);
                          }
        }
      }

      read_junk(fpinput,fpecho);
      for (i=0; i<Nwall_type; i++){
        if (Mix_type==0) {jmin=i; jmax=i+1;}
        else if (Mix_type==1) {jmin=0;jmax=Nwall_type;}
        for (j=jmin; j<jmax; j++){
	  if (Mix_type==1){ fscanf(fpinput,"%lf",&YukawaK_ww[i][j]);
	                     fprintf(fpecho,"%f  ",YukawaK_ww[i][j]);
                          }
          else            { fscanf(fpinput,"%lf",&YukawaK_w[i]);
	                     fprintf(fpecho,"%f  ",YukawaK_w[i]);
                          }
        }
      }
  }
  else{
      read_junk(fpinput,fpecho);  fprintf(fpecho,"EpsYukawa_ww n/a");
      read_junk(fpinput,fpecho);  fprintf(fpecho,"YukawaK_ww n/a");
  }

  /*********************************************************/
  /* Initialize and Read Wall-Fluid Interaction Parameters */
  /*********************************************************/

  if (Mix_type == 1){

     read_junk(fpinput,fpecho);
                       /* Sigma_wf */
     for (i=0; i<Ncomp; i++){
          for (j=0; j<Nwall_type; j++){
              fscanf(fpinput,"%lf",&Sigma_wf[i][j]);
               fprintf(fpecho,"%f  ",Sigma_wf[i][j]);
          }
     }
                       /* Eps_wf */
     read_junk(fpinput,fpecho);
     for (i=0; i<Ncomp; i++){
          for (j=0; j<Nwall_type; j++){
              fscanf(fpinput,"%lf",&Eps_wf[i][j]);
               fprintf(fpecho,"%f  ",Eps_wf[i][j]);
          }
     }
                       /* Cut_wf */
     read_junk(fpinput,fpecho);
     for (i=0; i<Ncomp; i++){
          for (j=0; j<Nwall_type; j++){
   	      fscanf(fpinput,"%lf",&Cut_wf[i][j]);
   	       fprintf(fpecho,"%f  ",Cut_wf[i][j]);
          }
     }
                       /* Yukawa Energy paramter */
     read_junk(fpinput,fpecho);
     for (i=0; i<Ncomp; i++){
          for (j=0; j<Nwall_type; j++){
   	      fscanf(fpinput,"%lf",&EpsYukawa_wf[i][j]);
   	       fprintf(fpecho,"%f  ",EpsYukawa_wf[i][j]);
          }
     }
                     /* Yukawa paramter */
     read_junk(fpinput,fpecho);
     for (i=0; i<Ncomp; i++){
         for (j=0; j<Nwall_type; j++){
 	      fscanf(fpinput,"%lf",&YukawaK_wf[i][j]);
 	       fprintf(fpecho,"%f  ",YukawaK_wf[i][j]);
         }
     }
  }
  else{
     read_junk(fpinput,fpecho);
      fprintf(fpecho,"\n  USING L_B Mixing Rules --- WALL-FLUID INTERACTIONS COMPUTED BY CODE\n");
      fprintf(fpecho,"........MANUAL INPUT DOES NOT APPLY \n");
      fprintf(fpecho,"skipping this parameter  ");
     for (i=0; i<4; i++){
          read_junk(fpinput,fpecho);
          fprintf(fpecho,"skipping this parameter  ");
     }
  }

  /************************************************/
  /* Initialize and Read Polymer Input Parameters */
  /************************************************/

  for (icomp=0;icomp<Ncomp;icomp++) Unk2Comp[icomp]=icomp;    /* this will be reset for polymers where we track segments */

  if (Type_poly != NONE){

    read_junk(fpinput,fpecho);
    fscanf(fpinput,"%d",&Npol_comp);
     fprintf(fpecho,"%d",Npol_comp);

    read_junk(fpinput,fpecho);
    for (i=0; i<Npol_comp; ++i){
	fscanf(fpinput,"%d",&Nblock[i]);
	 fprintf(fpecho,"%d",Nblock[i]);
	if (Nblock[i] > NBLOCK_MAX) {
	  printf("Error: Must increase NBLOCK_MAX");
	  exit(-1);
        }
    }

    read_junk(fpinput,fpecho);
    for (pol_number=0; pol_number<Npol_comp; ++pol_number){
	Nmer[pol_number] = 0;
	for (i=0; i<Nblock[pol_number]; ++i){
	  fscanf(fpinput,"%d", &block[pol_number][i]);
	   fprintf(fpecho,"%d  ",block[pol_number][i]);
          Nseg_per_block[pol_number][i]=block[pol_number][i];
	  Nmer[pol_number] += block[pol_number][i];
	}
    }

    Ntype_mer = 0;
    read_junk(fpinput,fpecho);
    seg_tot=0;
    for (i=0; i<NBLOCK_MAX; ++i){ 
        Nmer_t_total[i]=0; 
        Type_mer_to_Pol[i]=-1;
    }
    for (pol_number=0; pol_number<Npol_comp; ++pol_number){
         Poly_to_Ntype[pol_number]=0;
         seg = 0;
         for (i=0; i<NBLOCK_MAX; ++i)  Nmer_t[pol_number][i] = 0; 

	 for (i=0; i<Nblock[pol_number]; ++i){
	     fscanf(fpinput,"%d", &block_type[i]);
	      fprintf(fpecho,"%d  ",block_type[i]);
             SegType_per_block[pol_number][i]=block_type[i];
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

      graft_logical=FALSE;
      read_junk(fpinput,fpecho);
      for (pol_number=0; pol_number<Npol_comp; ++pol_number){
          fscanf(fpinput,"%d",&Grafted[pol_number]);
           fprintf(fpecho,"%d   ",Grafted[pol_number]);
          if (Grafted[pol_number]==TRUE) graft_logical=TRUE;
      }

      read_junk(fpinput,fpecho);
      for (pol_number=0; pol_number<Npol_comp; ++pol_number){
         if(graft_logical==TRUE) {
            fscanf(fpinput,"%d", &Graft_wall[pol_number]);
             fprintf(fpecho,"%d  ",Graft_wall[pol_number]);
         }
         else Graft_wall[pol_number]=-1;
      }

      read_junk(fpinput,fpecho);
      for (pol_number=0; pol_number<Npol_comp; ++pol_number){
          if(graft_logical==TRUE) {
              fscanf(fpinput,"%lf",&Rho_g[pol_number]);
               fprintf(fpecho,"%f   ",Rho_g[pol_number]);
          }
          else Rho_g[pol_number]=0.0;
      }

      read_junk(fpinput,fpecho);
      fscanf(fpinput,"%d", &Type_poly_arch);
      if (Type_poly_arch==POLY_ARCH_FILE && Type_poly != SCFT){
           fscanf(fpinput,"%s", poly_file_array);
            fprintf(fpecho,"%s  ",poly_file_array);
           Poly_file_name=poly_file_array;
      }

    /* now set up the chain architecture...either from a file or from known structure */
    setup_chain_architecture(poly_file_array,fpecho);

    if (Type_poly == CMS){  /* this bit only applies to the CMS functional */
       read_junk(fpinput,fpecho);
       fscanf(fpinput,"%d",&Ncr_files);
       fscanf(fpinput,"%lf", &Crfac);
       fscanf(fpinput,"%s", cr_file_array);
        fprintf(fpecho,"%d  %9.6f  %s ", Ncr_files,Crfac,cr_file_array);
       Cr_file=cr_file_array;

       if (Ncr_files>=2){ 
             fscanf(fpinput,"%s",cr_file2_array);
              fprintf(fpecho,"  %s ",cr_file2_array);
             Cr_file2=cr_file2_array;
       }
/*       if (Ncr_files>=3){ 
             fscanf(fpinput,"%s",Cr_file3);
              fprintf(fpecho,"  %s ",Cr_file3);
       }
       if (Ncr_files>=4){ 
             fscanf(fpinput,"%s",Cr_file4);
              fprintf(fpecho,"  %s ",Cr_file4);
       }*/
       read_junk(fpinput,fpecho);
/*        if ( fabs(Crfac+1.0)<1.e-8){*/
         for (i=0;i<Ncr_files-2;i++){
            fscanf(fpinput,"%lf", &Cr_break[i]);
             fprintf(fpecho,"%f  ",Cr_break[i]);
         }
        /*}
       else  fprintf(fpecho,"n/a: not doing automated linear interpolation");*/
       read_junk(fpinput,fpecho);
       for  (icomp=0; icomp<Ncomp; icomp++){
          for  (jcomp=0; jcomp<Ncomp; jcomp++){
             fscanf(fpinput,"%lf", &Cr_rad_hs[icomp][jcomp]);
              fprintf(fpecho,"%f  ",Cr_rad_hs[icomp][jcomp]);
          }
       }
    /* Note: the value of Cr_rad_hs may get reset to Cutoff_ff if Ipot_ff_n=2 
       see setup_polymer_cr in dft_main.c */
       }
       else{
          read_junk(fpinput,fpecho);
           fprintf(fpecho,"\n NO LIQUID STATE INPUT FOR NON-CMS RUN\n  n/a   ");
          for (i=0; i<2; i++) { read_junk(fpinput,fpecho);  fprintf(fpecho,"n/a   "); }       
       }
  }    
  else{
     read_junk(fpinput,fpecho);
      fprintf(fpecho,"\n POLYMER INPUT NOT RELEVENT FOR THIS RUN\n not read  ");
     for (i=0; i<10; i++) { read_junk(fpinput,fpecho);  fprintf(fpecho,"not read  - no polymers  "); }
  }
  /* end of polymer input */

  /***************************************************************/
  /* Initialize and Read Parameters for Semi-Permeable Membranes */
  /***************************************************************/

  if (Nwall_type >0) Lsemiperm = (int **) array_alloc (2, Nwall_type,Ncomp,sizeof(int));
  lzeros=FALSE;
  ltrues=FALSE;
  read_junk(fpinput,fpecho);
  for (i=0; i<Nwall_type; i++){
      for(j=0; j<Ncomp; j++){
        if (lzeros==FALSE  && ltrues==FALSE){
           fscanf(fpinput,"%d ", &Lsemiperm[i][j]);
            fprintf(fpecho,"%d  ",Lsemiperm[i][j]);
        }
        if (i==0 && j==0){
           if (Lsemiperm[i][j]==-2) lzeros=TRUE;
           else if (Lsemiperm[i][j]==-1) ltrues=TRUE;
        }
        if (lzeros==TRUE) Lsemiperm[i][j]=0;
        else if (ltrues==TRUE) Lsemiperm[i][j]=1;
      }
  }

  if (Nwall_type>0) Vext_membrane = (double **) array_alloc (2, Nwall_type,Ncomp,sizeof(double));
  read_junk(fpinput,fpecho);
  for (i=0; i<Nwall_type; i++){
      for(j=0; j<Ncomp; j++){
        if (lzeros){ Vext_membrane[i][j]=0.0; }
        else{
          fscanf(fpinput,"%lf", &Vext_membrane[i][j]);
           fprintf(fpecho,"%f  ",Vext_membrane[i][j]);
        }
      }
  }
  for (i=0; i<Nwall_type; i++){
        for(j=0; j<Ncomp; j++) if (Lsemiperm[i][j]==FALSE) Vext_membrane[i][j] = VEXT_MAX;
  }

  /**********************************************/
  /* Initialize and Read State Point Parameters */
  /**********************************************/

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d  %d %d", &Type_interface, &Grad_dim,&Lconstrain_interface);
   fprintf(fpecho,"%d  %d  %d",Type_interface,Grad_dim,Lconstrain_interface);

   read_junk(fpinput,fpecho);
   if (Type_poly==NONE || Ntype_mer == 1){
      for (icomp=0; icomp<Ncomp; ++icomp){
        fscanf(fpinput,"%lf", &Rho_b[icomp]);
         fprintf(fpecho,"%f  ",Rho_b[icomp]);
        if (Type_interface != UNIFORM_INTERFACE) Rho_b_LBB[icomp]=Rho_b[icomp];
      }
   }
   else{
      /*first scan in polymer densities */
      for (i=0; i<Npol_comp; i++) fscanf(fpinput,"%lf", &rho_tmp[i]);
      for (i=0; i<Ntype_mer; i++){
        Rho_b[i] = 0.;
        for (j=0; j<Npol_comp; j++)  Rho_b[i] += (double)Nmer_t[j][i]*rho_tmp[j]/(double)Nmer[j];

        if (Type_interface != UNIFORM_INTERFACE) Rho_b_LBB[i]=Rho_b[i];
         fprintf(fpecho,"%f  ",Rho_b[i]);
      }
  }

  read_junk(fpinput,fpecho);
  if (Type_interface != UNIFORM_INTERFACE){
       if (Type_poly==NONE || Ntype_mer == 1){
         for (icomp=0; icomp<Ncomp; ++icomp){
           fscanf(fpinput,"%lf", &Rho_b[icomp]);
            fprintf(fpecho,"%f  ",Rho_b[icomp]);
           if (Type_interface != UNIFORM_INTERFACE) Rho_b_RTF[icomp]=Rho_b[icomp];
         }
       }
       else{
         /*first scan in polymer densities */
         for (i=0; i<Npol_comp; i++){
             fscanf(fpinput,"%lf", &rho_tmp[i]);
         }
         for (i=0; i<Ntype_mer; i++){
           Rho_b[i] = 0.;
           for (j=0; j<Npol_comp; j++)  Rho_b[i] += (double)Nmer_t[j][i]*rho_tmp[j]/(double)Nmer[j];

           if (Type_interface != UNIFORM_INTERFACE) Rho_b_RTF[i]=Rho_b[i];
            fprintf(fpecho,"%f  ",Rho_b[i]);
         }
       }
  }
  else{  fprintf(fpecho,"n/a"); }

  read_junk(fpinput,fpecho);
  if (Ipot_ff_c == COULOMB && Type_interface != UNIFORM_INTERFACE) {
     fscanf(fpinput,"%lf", &Elec_pot_LBB);
      fprintf(fpecho,"%f  ",Elec_pot_LBB);
     fscanf(fpinput,"%lf", &Elec_pot_RTF);
      fprintf(fpecho,"%f  ",Elec_pot_RTF);
  }
  else{  fprintf(fpecho,"n/a"); }

  read_junk(fpinput,fpecho);
  if (Type_interface != UNIFORM_INTERFACE){
         fscanf(fpinput,"%lf", &X_const_mu);
          fprintf(fpecho,"%f  ",X_const_mu);
  }
  else{  fprintf(fpecho,"n/a"); }


  /**************************************************************/
  /* Initialize and Read Parameters Specific to Charged Systems */
  /**************************************************************/
  Ipot_wf_c = 0;
  if (Type_coul != NONE || Type_pairPot==PAIR_COULOMB_CS || Type_pairPot==PAIR_COULOMB) {
     read_junk(fpinput,fpecho);
     for (iwall_type=0; iwall_type <  Nwall_type; ++iwall_type){
	fscanf(fpinput,"%d", &Type_bc_elec[iwall_type]);
	 fprintf(fpecho,"%d  ",Type_bc_elec[iwall_type]);
        if (Type_bc_elec[iwall_type] != 0) Ipot_wf_c = COULOMB;
     }

     read_junk(fpinput,fpecho);
     fscanf(fpinput,"%d", &Nlocal_charge);
      fprintf(fpecho,"%d  ",Nlocal_charge);
     if (Nlocal_charge > 0) Ipot_wf_c=COULOMB;

     ncharge=0;
     if      (Nlocal_charge > 0)  ncharge = Nlocal_charge;
     else if (Nlocal_charge < 0)  ncharge = 2;
    
     read_junk(fpinput,fpecho);
     if (Nlocal_charge == 0)  {
           fprintf(fpecho,"n/a: Nlocal_charge=0");
          read_junk(fpinput,fpecho);
           fprintf(fpecho,"n/a: Nlocal_charge=0");
          read_junk(fpinput,fpecho);
           fprintf(fpecho,"n/a: Nlocal_charge=0");
     }
     else{
       Charge      = (double *) array_alloc (1, ncharge, sizeof(double));
       Charge_Diam = (double *) array_alloc (1, ncharge, sizeof(double));
       Charge_x = (double **) array_alloc (2, Ndim,ncharge,sizeof(double));
       for(i=0; i<ncharge; i++){
	  fscanf(fpinput,"%lf ", &Charge[i]);
	   fprintf(fpecho,"%f  ",Charge[i]);
       }
       read_junk(fpinput,fpecho);
       for(i=0; i<ncharge; i++){
	  fscanf(fpinput,"%lf ", &Charge_Diam[i]);
	   fprintf(fpecho,"%f  ",Charge_Diam[i]);
       }
	read_junk(fpinput,fpecho);
	for(i=0; i<ncharge; i++){
	  for (idim=0; idim <  Ndim; ++idim){
	    fscanf(fpinput,"%lf", &Charge_x[idim][i]);
	     fprintf(fpecho,"%f  ",Charge_x[idim][i]);
	  }
	}
    }
    logical=FALSE;
    for (iwall_type=0; iwall_type<Nwall_type; iwall_type++)
      if (Type_bc_elec[iwall_type] == ATOMIC_CHARGE) logical=TRUE;

    read_junk(fpinput,fpecho);
    if (ncharge !=0 || logical){
      fscanf(fpinput,"%d  %d", &Charge_type_atoms, &Charge_type_local);
       fprintf(fpecho,"%d  %d ",Charge_type_atoms,Charge_type_local);
    }
    else  fprintf(fpecho,"n/a:  ncharge=0 and walls are not atoms");
  }
  else{
      read_junk(fpinput,fpecho);
       fprintf(fpecho,"CHARGED SURFACE INPUT NOT RELEVENT FOR THIS RUN\n");
       fprintf(fpecho,"n/a");
      for (i=0; i<5; i++){ read_junk(fpinput,fpecho);  fprintf(fpecho,"n/a"); }
  }
  
   /* Dielectric constant choices */

  if (Ipot_ff_c == COULOMB ) {
     read_junk(fpinput,fpecho);
     fscanf(fpinput,"%d  %lf %lf %lf ",&Type_dielec, &Sigma_Angstroms_plasma, &Temp_K_plasma, &DielecConst_plasma);
      fprintf(fpecho,"%d  %f %f %f",Type_dielec,Sigma_Angstroms_plasma,Temp_K_plasma,DielecConst_plasma);

     read_junk(fpinput,fpecho);
     if (Type_dielec != DIELEC_WF_PORE){
       fscanf(fpinput,"%lf",&Dielec_bulk);
        fprintf(fpecho,"%f  ",Dielec_bulk);
     }
     else {
       fscanf(fpinput,"%lf %lf %lf",&Dielec_bulk, &Dielec_pore, &Dielec_X);
        fprintf(fpecho,"%f  %f  %f ",Dielec_bulk,Dielec_pore,Dielec_X);
     }
    
     read_junk(fpinput,fpecho);
     if (Nwall_type > 0) {
       Dielec_wall = (double *) array_alloc (1, Nwall_type,sizeof(double));
       for (iwall_type=0; iwall_type<Nwall_type; iwall_type++) {
	  fscanf(fpinput,"%lf", &Dielec_wall[iwall_type]);
	   fprintf(fpecho,"%f  ",Dielec_wall[iwall_type]);
       }
     }
     else  fprintf(fpecho,"n/a");
  }
  else{
      read_junk(fpinput,fpecho);
       fprintf(fpecho,"CHARGED FLUID (DIELECTRIC) INPUT NOT RELEVENT FOR THIS RUN\n");
       fprintf(fpecho,"n/a");
      for (i=0; i<2; i++){
	read_junk(fpinput,fpecho);
	 fprintf(fpecho,"n/a");
      }
  }
  
  /**************************************************************/
  /* Initialize and Read Parameters Specific to Diffusing Systems */
  /**************************************************************/

  L1D_bc = FALSE;
  Linear_transport=0;   /* only want generalized Fick's Law */
                        /* can do J=-Dgrad mu, but doesn't work well */

  if (Type_interface==DIFFUSIVE_INTERFACE) {

    read_junk(fpinput,fpecho);
    for (icomp=0; icomp<Ncomp; ++icomp){
	fscanf(fpinput,"%lf", &D_coef[icomp]);
	 fprintf(fpecho,"%f  ",D_coef[icomp]);
    }

    read_junk(fpinput,fpecho);
    fscanf(fpinput,"%lf ", &Velocity);
     fprintf(fpecho,"%f  ",Velocity);

    read_junk(fpinput,fpecho);
    fscanf(fpinput,"%d  %d", &Geom_flag, &Nseg_IC);
     fprintf(fpecho,"%d  %d  ",Geom_flag,Nseg_IC);

    Pore_rad_L_IC = (double *) array_alloc (1, Nseg_IC,sizeof(double));
    Pore_rad_R_IC = (double *) array_alloc (1, Nseg_IC,sizeof(double));
    Lseg_IC       = (double *) array_alloc (1, Nseg_IC,sizeof(double));
    read_junk(fpinput,fpecho);
      for (i=0; i<Nseg_IC; i++){
	fscanf(fpinput,"%lf  %lf  %lf", &Pore_rad_L_IC[i], &Pore_rad_R_IC[i], &Lseg_IC[i]);
	 fprintf(fpecho,"%f  %f  %f  ", Pore_rad_L_IC[i],Pore_rad_R_IC[i],Lseg_IC[i]);
        if (Length_ref > 0.0){
            Pore_rad_L_IC[i]/=Length_ref;
            Pore_rad_R_IC[i]/=Length_ref;
            Lseg_IC[i]/=Length_ref;
        }
    }

  }
  else{
      read_junk(fpinput,fpecho);
       fprintf(fpecho,"TRANSPORT INPUT NOT RELEVENT FOR THIS RUN\n n/a ");
      for (i=0; i<3; i++) { read_junk(fpinput,fpecho);  fprintf(fpecho,"n/a"); }
  }

  /******************************************************************/
  /* Initialize and Read Settings for Run control and Initial Guess */
  /******************************************************************/

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d %d",&Iguess,&Iguess_fields);
   fprintf(fpecho,"%d  %d ",Iguess,Iguess_fields);

  if (Iguess != STEP_PROFILE) Nsteps=1;
  if (Iguess==STEP_PROFILE || (Iguess>=CHOP_RHO && Iguess<= CHOP_RHO_STEP)){
      read_junk(fpinput,fpecho);
      fscanf(fpinput,"%d",&Nsteps);
       fprintf(fpecho,"%d   ",Nsteps);

      read_junk(fpinput,fpecho);
      for (i=0; i<Nsteps; ++i){
	fscanf(fpinput,"%d", &Orientation_step[i]);
	 fprintf(fpecho,"%d  ",Orientation_step[i]);
      }

      read_junk(fpinput,fpecho);
      for (i=0; i<Nsteps; ++i){
	fscanf(fpinput,"%lf", &Xstart_step[i]);
	 fprintf(fpecho,"%f  ",Xstart_step[i]);
        Xstart_step[i]-=1.0e-4; /* prevent roundoff errors later*/
      }

      read_junk(fpinput,fpecho);
      for (i=0; i<Nsteps; ++i){
	fscanf(fpinput,"%lf", &Xend_step[i]);
	 fprintf(fpecho,"%f  ",Xend_step[i]);
        Xend_step[i]+=1.e-4; /* prevent roundoff errors later*/
      }

      read_junk(fpinput,fpecho);
      for (icomp=0; icomp<Ncomp; ++icomp){
          for (i=0; i<Nsteps; ++i){
   	    fscanf(fpinput,"%lf", &Rho_step[icomp][i]);
	     fprintf(fpecho,"%f  ",Rho_step[icomp][i]);
          }
      }
  }
  else{ for (i=0;i<5;i++) read_junk(fpinput,fpecho); }

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d ",&Restart);
  if (Restart == RESTART_FEWERCOMP) fscanf(fpinput," %d",&Nmissing_densities);
  else Nmissing_densities=0;
   fprintf(fpecho,"%d  %d",Restart,Nmissing_densities);

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d ",&Restart_Vext);
   fprintf(fpecho,"%d ",Restart_Vext);
  if (Restart_Vext != READ_VEXT_FALSE){
        fscanf(fpinput,"%s", vext_file_array);
         fprintf(fpecho,"  %s ", vext_file_array);
           Vext_filename=vext_file_array;
        if (Restart_Vext == READ_VEXT_SUMTWO || Restart_Vext == READ_VEXT_STATIC){
           fscanf(fpinput,"%s", vext_file2_array);
            fprintf(fpecho,"  %s ", vext_file2_array);
           Vext_filename2=vext_file2_array;
        }
        else Vext_filename2=NULL;
  }
  else{
     Vext_filename=NULL;
     Vext_filename2=NULL;
  }

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%lf",&Rho_max);
   fprintf(fpecho,"%f  ",Rho_max);

  /*********************************************/
  /* Initialize and Read Output Format Setting */
  /*********************************************/

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d  %d  %d %d ",&Lper_area, &Lcount_reflect, &Lprint_gofr, &Lprint_pmf);
   fprintf(fpecho,"%d  %d  %d  %d",Lper_area,Lcount_reflect,Lprint_gofr,Lprint_pmf);

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d",&Print_rho_type);
   fprintf(fpecho,"%d",Print_rho_type);
  if (Print_rho_type==0 && Restart != 0){
      printf("WARNING: Print_rho_type is being set to 1 so that restart files will not be overwritten\n");
      Print_rho_type=1;
  }
  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d %d",&Print_rho_switch,&Print_mesh_switch);
   fprintf(fpecho,"%d  %d",Print_rho_switch,Print_mesh_switch);
  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d",&Iwrite);
   fprintf(fpecho,"%d",Iwrite);

  if (Iwrite==MINIMAL || Iwrite==DENSITIES){
    Iwrite_screen=SCREEN_BASIC;
    Iwrite_files=FILES_BASIC;
  }
  else if (Iwrite==EXTENDED){
    Iwrite_screen=SCREEN_BASIC;
    Iwrite_files=FILES_EXTENDED;
  }
  else if (Iwrite==VERBOSE){
    Iwrite_screen=SCREEN_VERBOSE;
    Iwrite_files=FILES_DEBUG;
  }
  else if (Iwrite==NO_SCREEN){
    Iwrite_screen=SCREEN_NONE;
    Iwrite_files=FILES_BASIC;
  }
  else if (Iwrite==VERBOSE_MATRIX){
    Iwrite_screen=SCREEN_BASIC;
    Iwrite_files=FILES_DEBUG_MATRIX;
  }

  if (fabs(charge_sum) > 1.e-8 && Iwrite_screen != SCREEN_NONE && Iwrite_screen != SCREEN_ERRORS_ONLY) printf("\n TOTAL CHARGE IN dft_surfaces.dat = %9.6f\n",charge_sum);

  /**************************************************/
  /* Numerical Methods ... Mesh/Jacobian Coarsening */
  /**************************************************/

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d",&Nzone);
   fprintf(fpecho,"%d",Nzone);

  read_junk(fpinput,fpecho);
  if (Nzone != 1) {
    for (izone =0; izone<Nzone-1; izone++){
        fscanf(fpinput,"%lf",&Rmax_zone[izone]);
         fprintf(fpecho,"%f  ",Rmax_zone[izone]);
    }
  }

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d",&Mesh_coarsening);
   fprintf(fpecho,"%d",Mesh_coarsening);

  if (Mesh_coarsening==2 && Lauto_size){
      if(Iwrite_screen != SCREEN_NONE && Iwrite_screen != SCREEN_ERRORS_ONLY) printf("AUTO SIZE FOR BOX :: [Lx,Ly,Lz]=");
      for (idim=0; idim<Ndim; idim++){
         Size_x[idim] = maxpos[idim]-minpos[idim] + 2.0*(Rmax_zone[0]+Sigma_ff[0][0]);
         if(Iwrite_screen != SCREEN_NONE && Iwrite_screen != SCREEN_ERRORS_ONLY) printf("  %9.6f\n",Size_x[idim]);
      }
  }

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d %lf",&Coarser_jac,&Jac_grid);
   fprintf(fpecho,"%d %f ",Coarser_jac,Jac_grid);

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d %lf",&Lcut_jac,&Jac_threshold);
   fprintf(fpecho,"%d %f ",Lcut_jac,Jac_threshold);

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d  %lf", &L1D_bc,&X_1D_bc );
   fprintf(fpecho,"%d  %f  ",L1D_bc,X_1D_bc);
  
  if (L1D_bc==TRUE) Dim_1Dbc=Grad_dim; /* old code did this - GUI has new parameter */
  
  /*****************************************************/
  /* Numerical Methods ... Nonlinear Solver Parameters */
  /*****************************************************/

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d %d %d %d %d", &NL_Solver, &Max_NL_iter, &Physics_scaling, &ATTInA22Block, &Analyt_WJDC_Jac);
   fprintf(fpecho,"NL_Solver=%d %d %d %d %d",NL_Solver,Max_NL_iter,Physics_scaling,ATTInA22Block,Analyt_WJDC_Jac);
  read_junk(fpinput,fpecho);
  if (Physics_scaling != FALSE  && (Type_poly == WJDC || Type_poly==WJDC2 || Type_poly==WJDC3)){
      for (ipol=0;ipol<Npol_comp;ipol++){
            for (icomp=0;icomp<Ncomp;icomp++) {
                  fscanf(fpinput,"%lf ", &Scale_fac_WJDC[ipol][icomp]);
                   fprintf(fpecho,"%f ", Scale_fac_WJDC[ipol][icomp]);
            }
      }
  }
  else { if(Iwrite_screen == SCREEN_VERBOSE) printf("n/a - scaling parameters not relevant"); }

  if (NL_Solver==PICARD_BUILT_IN && Iguess_fields !=CALC_ALL_FIELDS){
     if(Iwrite_screen != SCREEN_NONE) printf("Picard solver indicated so Iguess_fields is reset to %d\n",CALC_ALL_FIELDS);
     Iguess_fields=CALC_ALL_FIELDS;
  }

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%lg", &NL_rel_tol);
   fprintf(fpecho,"%lg  ",NL_rel_tol);
  fscanf(fpinput,"%lg", &NL_abs_tol);
   fprintf(fpecho,"%lg  ",NL_abs_tol);
  fscanf(fpinput,"%lg", &NL_update_scalingParam);
   fprintf(fpecho,"%lg  ",NL_update_scalingParam);

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%lg", &NL_rel_tol_picard);
   fprintf(fpecho,"%lg  ",NL_rel_tol_picard);
  fscanf(fpinput,"%lg", &NL_abs_tol_picard);
   fprintf(fpecho,"%lg  ",NL_abs_tol_picard);

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d", &Load_Bal_Flag);
   fprintf(fpecho,"%d  ", Load_Bal_Flag);

  /*****************************************************/
  /* Numerical Methods ... Linear Solver Parameters    */
  /*****************************************************/

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d ", &L_Schur);
  fscanf(fpinput,"%d ", &Az_solver);
  if (Az_solver == 0) fscanf(fpinput,"%d ", &Az_kspace);
  else Az_kspace=-1;
   fprintf(fpecho,"L_Schur=%d  %d  %d ", L_Schur, Az_solver, Az_kspace);

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d", &Az_scaling);
   fprintf(fpecho,"%d  ",Az_scaling);

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d  %lf", &Az_preconditioner, &Az_ilut_fill_param);
   fprintf(fpecho,"%d  %f ",Az_preconditioner, Az_ilut_fill_param);

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d  %lg", &Max_gmres_iter,&Az_tolerance);
   fprintf(fpecho,"%d  %g  ",Max_gmres_iter,Az_tolerance);

  /********************************************************/
  /* Initialize and Read Parameters for Mesh Continuation */
  /********************************************************/

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d",&Nruns);
   fprintf(fpecho,"Nruns=%d  ",Nruns);

  if (Nruns > 0){

    read_junk(fpinput,fpecho);
    for (i=0; i < Ndim; ++i){
      fscanf(fpinput,"%lf", &Del_1[i]);
       fprintf(fpecho,"%f  ",Del_1[i]);
    }

    read_junk(fpinput,fpecho);
    fscanf(fpinput,"%d",&Plane_new_nodes);
    fscanf(fpinput,"%d",&Pos_new_nodes);
     fprintf(fpecho,"%d  %d ",Plane_new_nodes,Pos_new_nodes);
  }
  else{ read_junk(fpinput,fpecho); read_junk(fpinput,fpecho);}

  /********************************************************/
  /* Initialize and Read Parameters for LOCA Continuation */
  /********************************************************/

#ifdef USE_LOCA
  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d", &itmp);
   fprintf(fpecho,"%d ",itmp);
  Loca.method = itmp;
/*   fprintf(fpecho,"%d  ",Loca.method);*/

  if (Loca.method == 4) Lbinodal = 1;
  else Lbinodal=0;

  if (Lbinodal && Restart==0){
   printf("ERROR: Can only do binodal calculations with a restart\n") ; exit (-1);
  }

  for (i=0;i<NCONT_MAX;i++){ for (j=0;j<2;j++) Cont_ID[i][j]=-1; } 
  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d %d ", &itmp, &NID_Cont);
  for (i=0;i<NID_Cont;i++) fscanf(fpinput,"%d ", &Cont_ID[0][i]);
  Loca.cont_type1 = itmp;

   fprintf(fpecho,"%d  %d  ",Loca.cont_type1,NID_Cont);
  for (i=0;i<NID_Cont;i++)  fprintf(fpecho,"%d  ",Cont_ID[0][i]);

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%lf", &dtmp);

  if (Loca.method !=-1) {
    if ( Loca.cont_type1 == CONT_RHO_I || 
              Loca.cont_type1 == CONT_LOG_RHO_I){
           if (Type_poly == CMS){
               printf("WARNING: density continuation may be incorrect with polymer code\n");
               printf("...... you really need a new c(r) for the hard chain part\n");
               printf("...... for each new density.\n");
           }
    }
    if (Mix_type == 0 && (Loca.cont_type1 == CONT_EPSWF_IJ ||  
                          Loca.cont_type1 == CONT_EPSWF_ALL)){
        printf("ERROR: Can't do continuation in Eps_wf when the Mix_type is 0\n");
        exit(-1);
    } 
  }

  Loca.step_size = dtmp;
   fprintf(fpecho,"%f  ",Loca.step_size);

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d", &itmp);
  Loca.num_steps = itmp;

  fscanf(fpinput,"%lf", &dtmp);
  Loca.aggr = dtmp;
   fprintf(fpecho,"%d %f ",Loca.num_steps, Loca.aggr);

  read_junk(fpinput,fpecho);
  fscanf(fpinput,"%d %d ", &itmp, &NID_Cont);
  for (i=0;i<NID_Cont;i++) fscanf(fpinput,"%d ", &Cont_ID[1][i]);
  Loca.cont_type2 = itmp;

   fprintf(fpecho,"%d  %d  ",Loca.cont_type2,NID_Cont);
  for (i=0;i<NID_Cont;i++)  fprintf(fpecho,"%d  ",Cont_ID[1][i]);

  if (!LBulk && (((Loca.method != -1 && Loca.cont_type1 == CONT_BETAMU_I) ||
       (Loca.method==4 && Loca.cont_type2 == CONT_BETAMU_I)) 
      || ((Loca.method != -1 && Loca.cont_type1 == CONT_BETAMU_I_NEW) ||
          (Loca.method==4 && Loca.cont_type2 == CONT_BETAMU_I_NEW))) 
     ){
       /*printf("for continuation in chemical potential LBulk must be TRUE=%d .... resetting LBulk\n",Loca.cont_type1,1);*/
       LBulk=TRUE;
    }

  /* checks on LBulk */
  /* first check that bulk boundaries are NOT used if LBulk=TRUE and chemical potentials will be varied */ 
  if (LBulk && Loca.method != -1 &&
     (Loca.cont_type1 == CONT_BETAMU_I) ){
     for (idim=0;idim<Ndim;idim++){
         for (i=0;i<2;i++){
            if (Type_bc[idim][i]==IN_BULK){ 
               printf("Bulk boundary detected while LBulk=TRUE and continuation in chemical potential requested.\n");
               printf("This will not work because Rho_b is not updated during a continuation run. \n");
               printf("The boundary condition will be reset to REFLECT (2). \n");
               Type_bc[idim][i]=LAST_NODE;
               Type_bc[idim][i]=REFLECT;
            }
         }
     } 
  }

#else
  Loca.method = -1;
   fprintf(fpecho,"Loca library not compiled in: %d  ",Loca.method);
#endif

  error_check();

  return;
}

/****************************************************************************/

/* read_junk :  This reads and then prints everything (i.e. the junk)
                between the last real input and the next.  It prints
                this junk to a file so that we will have saved the
                input file for the run of interest.                 */
void read_junk(FILE *fpinput, FILE *fpecho)
{
   int c;
   fprintf(fpecho,"   junk:");
   while ((c=getc(fpinput)) != EOF && c !='@') putc(c,fpecho);
   fprintf(fpecho,"   data:");

   return;
}
/****************************************************************************/
/* error_check :  This function checks for error or inconsistencies in
                  the input parameters.  If either are found, the program
                  is stopped and error messages are sent to dft_errors.dat */

void error_check(void)
{
  double charge_b;
  int i;
  char *yo="error_check";

  if (Ndim> NDIM_MAX) {
    if(Iwrite_screen!=SCREEN_NONE) printf("%s: ERROR: Ndim (%d) larger than NDIM_MAX ", yo, Ndim);
    if(Iwrite_screen!=SCREEN_NONE) printf("(%d).\n\t\t Recompile with bigger NDIM_MAX\n", NDIM_MAX);
    exit(-1);
  }
  if (Nwall > NWALL_MAX){
     if(Iwrite_screen!=SCREEN_NONE) printf ("\n Nwall (%d) too large ... increase NWALL_MAX (%d)\n",Nwall,NWALL_MAX);
     exit(-1);
  }
  if (Nwall_type > NWALL_MAX_TYPE){
     if(Iwrite_screen!=SCREEN_NONE) printf ("\n Nwall_type (%d) too large ... increase NWALL_MAX_TYPE (%d)\n",Nwall,NWALL_MAX_TYPE);
     exit(-1);
  }
  if (Ncomp> NCOMP_MAX) {
    if(Iwrite_screen!=SCREEN_NONE) printf("%s: ERROR: Ncomp (%d) larger than NCOMP_MAX ", yo, Ncomp);
    if(Iwrite_screen!=SCREEN_NONE) printf("(%d).\n\t\t  Recompile with bigger NCOMP_MAX\n", NCOMP_MAX);
    exit(-1);
  }


  if (Ipot_ff_n > 2){
     if(Iwrite_screen!=SCREEN_NONE) printf ("\nSorry, your choice for the fluid-fluid interaction is not available\n");
     if(Iwrite_screen!=SCREEN_NONE) printf ("Ipot_ff_n: %d\n", Ipot_ff_n);
     exit (-1);
  }

  for (i=0;i<Nwall_type;i++){
     if (Ipot_wf_n[i] > 7 ){
        if(Iwrite_screen!=SCREEN_NONE) printf ("\nSorry, your choice for the wall-fluid interaction is not available\n");
        if(Iwrite_screen!=SCREEN_NONE) printf ("Ipot_wf_n[%d]: %d\n", i,Ipot_wf_n[i]);
        exit (-1);
     }
  }

  if ((Type_func > 3) || (Type_func < -1)){
     if(Iwrite_screen!=SCREEN_NONE) printf ("\nSorry, your choice for the hs functional is not yet available\n");
     if(Iwrite_screen!=SCREEN_NONE) printf ("Type_func: %d\n", Type_func);
     if(Iwrite_screen!=SCREEN_NONE) printf ("Reset Type_func to 0 for the Rosenfeld functional\n");
     exit (-1);
  }

  if (Type_attr > 0 && Type_attr < -1){
     if(Iwrite_screen!=SCREEN_NONE) printf ("\nSorry, your choice for the type of attractions is not yet available\n");
     if(Iwrite_screen!=SCREEN_NONE) printf ("Type_attr: %d\n", Type_attr);
     if(Iwrite_screen!=SCREEN_NONE) printf ("Reset Type_attr to 0 for the strict mean field approximation\n");
     exit (-1);
  }

  if (Type_coul > 4 && Type_coul < -1){
     if(Iwrite_screen!=SCREEN_NONE) printf ("\nSorry, your choice for the Coulomb functionals\n");
     if(Iwrite_screen!=SCREEN_NONE) printf ("Type_coul: %d is not available \n", Type_coul);
     exit (-1);
  }

  if (Type_poly > WJDC3 || Type_poly < NONE){
     if(Iwrite_screen!=SCREEN_NONE) printf ("\nSorry, your choice for the type of polymer functional\n");
     if(Iwrite_screen!=SCREEN_NONE) printf ("Type_poly: %d is not available\n", Type_poly);
     exit (-1);
  }

  charge_b = 0.0;
  if (Ipot_ff_c !=0) {
     for (i=0; i<Ncomp; i++) {
             charge_b += Rho_b[i]*Charge_f[i];
     }
     if (charge_b > 1.0e-10) {
        if(Iwrite_screen!=SCREEN_NONE) printf ("ERROR:input file: the fluid is not neutral\n");
        if(Iwrite_screen!=SCREEN_NONE) printf ("charge_b: %lf\n",charge_b);
        exit (-1);
     }
  }

  if (Load_Bal_Flag < 0 || Load_Bal_Flag > 5) {
      if(Iwrite_screen!=SCREEN_NONE) printf("\nERROR: the Load_Bal_Flag (%d) must 0-5 \n", Load_Bal_Flag);
      exit (-1);
  }
  if (Num_Proc == 1) Load_Bal_Flag = 0;

  if (Nzone > NZONE_MAX || Nzone <1 ){
     if(Iwrite_screen!=SCREEN_NONE) printf("\nERROR: Nzone out of range: Minimum val=1; Maximum val=%d: Current val=%d\n",NZONE_MAX,Nzone);
     exit(-1);
  }

}
