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
 *  FILE: dft_guess_restart.c
 *
 *  This file contains routines for performing a restart from old files
 *
 */

#include "dft_guess_restart.h"
 
/************************************************************/
void guess_restart_from_files(int start_no_info,int iguess,double **xInBox)
{
  char filename[20];
  double *x_new,fac;
  int iunk,i;

  x_new = (double *) array_alloc(1, Nnodes*Nunk_per_node, sizeof(double));

  if (Proc == 0) {  /* Proc 0 reads in the data file */

     if ( Imain_loop == 0){

         /* START FROM AN OLD FILE - otherwise if on second or greater step in a 
                                    continuation run, all of the _old variables
                                    were set in collect_xold (dft_output.c) */

         if (Lbinodal && iguess==BINODAL_FLAG) sprintf(filename,"dft_dens2.dat");
         else                                  sprintf(filename,"dft_dens.dat");

         Nodes_old = find_length_of_file(filename);

                     /* Modify Nodes_old for the special case where we will read in a 1D file, but set up an initial
                         guess for a 2D or 3D system */
         if (Restart==RESTART_1DTOND)  Nodes_old=Nnodes;

         if (Lbinodal && iguess==BINODAL_FLAG){
                  X2_old = (double *) array_alloc(1, Nodes_old*Nunk_per_node, sizeof(double));
                  for (i=0;i<Nodes_old*Nunk_per_node;i++) X2_old[i]=0.0;
         }
         else     {X_old = (double *) array_alloc(1, Nodes_old*Nunk_per_node, sizeof(double));
                  for (i=0;i<Nodes_old*Nunk_per_node;i++) X_old[i]=0.0;}

         read_in_a_file(iguess,filename); /* Get X_old from the file! */
     }
     else{ 
                      /* Here we are on n>first continuation step.  Have all necessary fields */
         for (i=0;i<NEQ_TYPE;i++) Restart_field[i]=TRUE;
     }

printf("Imain_loop=%d\n  Nodes_old=%d  Nodes=%d\n",Imain_loop,Nodes_old,Nnodes);
     if (Nodes_old != Nnodes) {     /* Profile must be modified in some way.  Number of nodes in file does
                                       not match number of nodes in the current problem */

          /* fix up to allow for some mixing of a bulk solution with a previously
                               converged solution --this has been disabled with fac=1*/
         fac=1.0;
         shift_the_profile(x_new,fac);
     }
     else{
         for (iunk=0; iunk<Nunknowns; iunk++){
              if (iguess==Iguess1)                       x_new[iunk] = X_old[iunk];
              else if (iguess==BINODAL_FLAG && Lbinodal) x_new[iunk] = X2_old[iunk];
         }
     }
     if (iguess==Iguess1) safe_free((void *) &X_old);
     else if (iguess==BINODAL_FLAG && Lbinodal) safe_free((void *) &X2_old);

  }

  MPI_Bcast (Restart_field,NEQ_TYPE,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast (&start_no_info,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Restart_field[DENSITY]==FALSE) {
      printf("can't do automatic restart without density fields yet\n");
      exit(-1);
  }
  communicate_profile(x_new,xInBox);
  check_zero_densities(xInBox);

  safe_free((void *) &x_new);
  return;
}
/**********************************************************/
/*find_length_of_file: here just establish a file length.*/
int find_length_of_file(char *filename)
{
  int c,nodes_old=0;
  FILE *fp;
  if (Proc==0) {
    if( (fp=fopen(filename,"r")) == NULL){
      printf("Can't open file %s\n", filename);
      exit(1);
    }

    nodes_old=0;
    while ( (c=getc(fp)) != EOF){
      if (c == '\n') {
	nodes_old++;
	c=getc(fp);/*so that a following blank line won't increment nodes_old*/
      }
    }
    fclose(fp);
  }
  return(nodes_old);
}
/**************************************************************/
/*read_in_a_file: In this routine, we read in a file containing
                  an old solution and use it for our initial guess */
void read_in_a_file(int iguess,char *filename)
{
  int c;
  int i,iunk,junk,idim, inode,itype_mer,ipol,iseg,index,dim_tmp,iunk_file;
  int ijk_old[3],ijk_old_max[3],open_now,ndim_max,node_start;
  int unk_in_file, unk_start_in_file[NEQ_TYPE],header,eq_type;
  int unk_to_eq_in_file[3*NCOMP_MAX+NMER_MAX+NMER_MAX*NMER_MAX+13];
  int convert_to_comp_densities, convert_to_seg_densities,jseg,icomp;
  double pos_old,pos_old0[3],tmp,x_tmp;
  char filename2[20];
  char unk_char[20];
  FILE *fp5=NULL,*fp6=NULL;

  convert_to_seg_densities=FALSE; 
  convert_to_comp_densities=FALSE; 


                    /* open the dft_dens.dat file */
   if( (fp5=fopen(filename,"r")) == NULL){
     printf("Can't open file %s\n", filename);
     exit(1);
   }

                   /* identify which unknowns are found in this file and in what order */
    header=0;
    iunk=0;
    unk_in_file=0;
    for (eq_type=0;eq_type<NEQ_TYPE;eq_type++){

       fgets(unk_char,20,fp5);
       if (strncmp(unk_char,"DENSITY",5)==0) {
             Restart_field[DENSITY]=TRUE;
             header++;
             unk_in_file+=Ncomp;
             unk_start_in_file[DENSITY]=iunk;
             for (i=0;i<Ncomp;i++) unk_to_eq_in_file[iunk++]=DENSITY;
             if (Lseg_densities) convert_to_seg_densities=TRUE;
             else                convert_to_seg_densities=FALSE;
       }
       if (strncmp(unk_char,"DENSSEG",5)==0) {
             Restart_field[DENSITY]=TRUE;
             header++;
             unk_in_file+=Nseg_tot;
             unk_start_in_file[DENSITY]=iunk;
             for (i=0;i<Nseg_tot;i++) unk_to_eq_in_file[iunk++]=DENSITY;
             if (Lseg_densities) convert_to_comp_densities=FALSE;
             else                convert_to_comp_densities=TRUE;
       }
       if (strncmp(unk_char,"MFEQ",5)==0) {
             Restart_field[MF_EQ]=TRUE;
             header++;
             unk_in_file+=Ncomp;
             unk_start_in_file[MF_EQ]=iunk;
             for (i=0;i<Ncomp;i++) unk_to_eq_in_file[iunk++]=MF_EQ;
       }
       else if (strncmp(unk_char,"POISSON",5)==0){
             Restart_field[POISSON]=TRUE;
             header++;
             unk_in_file+=1;
             unk_start_in_file[POISSON]=iunk;
             unk_to_eq_in_file[iunk++]=POISSON;
       }
       else if (strncmp(unk_char,"HSRHOBAR",5)==0){
             Restart_field[HSRHOBAR]=TRUE;
             header++;
             unk_in_file+=Nrho_bar;
             unk_start_in_file[HSRHOBAR]=iunk;
             for (i=0;i<Nrho_bar;i++) unk_to_eq_in_file[iunk++]=HSRHOBAR;
       }
       else if (strncmp(unk_char,"CMSFIELD",5)==0){
             Restart_field[CMS_FIELD]=TRUE;
             header++;
             unk_in_file+=Ncomp;
             unk_start_in_file[CMS_FIELD]=iunk;
             for (i=0;i<Ncomp;i++) unk_to_eq_in_file[iunk++]=CMS_FIELD;
       }
       else if (strncmp(unk_char,"WJDCFIELD",5)==0){
             Restart_field[WJDC_FIELD]=TRUE;
             header++;
             unk_in_file+=Ncomp;
             unk_start_in_file[WJDC_FIELD]=iunk;
             for (i=0;i<Ncomp;i++) unk_to_eq_in_file[iunk++]=WJDC_FIELD;
       }
       else if (strncmp(unk_char,"CAVWTC",5)==0){
             Restart_field[CAVWTC]=TRUE;
             header++;
             unk_in_file+=2;
             unk_start_in_file[CAVWTC]=iunk;
             for (i=0;i<2;i++) unk_to_eq_in_file[iunk++]=CAVWTC;
       }
       else if (strncmp(unk_char,"BONDWTC",5)==0){
             Restart_field[BONDWTC]=TRUE;
             header++;
             unk_in_file+=Nbonds;
             unk_start_in_file[BONDWTC]=iunk;
             for (i=0;i<Nbonds;i++) unk_to_eq_in_file[iunk++]=BONDWTC;
       }
       else if (strncmp(unk_char,"CHEMPOT",5)==0){
             Restart_field[DIFFUSION]=TRUE;
             header++;
             if (Lseg_densities) unk_in_file+=Nseg_tot;
             else              unk_in_file+=Ncomp;
             unk_start_in_file[DIFFUSION]=iunk;
             if (Lseg_densities) for (i=0;i<Nseg_tot;i++) unk_to_eq_in_file[iunk++]=DIFFUSION;
             else                for (i=0;i<Ncomp;i++) unk_to_eq_in_file[iunk++]=DIFFUSION;
       }
    }
    if (Iwrite != NO_SCREEN) printf("Number of unknowns in the file=%d\n",unk_in_file);

    if (Lsteady_state && Restart_field[DIFFUSION]==FALSE) 
         if (Proc==0 && Iwrite != NO_SCREEN)
           printf("there is no chemical potential data in the restart file\n");
    if (Type_coul != NONE && Restart_field[POISSON]==FALSE)
           printf("there is no electrostatic potential data in the restart file\n");
    if ((Type_poly == CMS || Type_poly==CMS_SCFT) && Restart_field[CMS_FIELD]==FALSE)
           printf("there is no CMS field data in the restart file\n");
    if ((Type_poly == WJDC) && Restart_field[WJDC_FIELD]==FALSE)
           printf("there is no WJDC field data in the restart file\n");
    if (L_HSperturbation && Restart_field[HSRHOBAR]==FALSE)
           printf("there is no Rosenfeld nonlocal density data in the restart file\n");
    if (Type_attr==MF_VARIABLE && Restart_field[MF_EQ]==FALSE)
           printf("there is no mean field attractive variable data in the restart file \n");
    if (Restart_field[DENSITY]==FALSE)
           printf("there is no density data in the restart file\n");

    fclose(fp5);
    if (Restart != RESTART_1DTOND) Nodes_old-=header;
    printf("skipping %d lines in the dft_dens.dat file\n",header);

  /* read positions from file find Nodes_x_old[idim] */
  /* read the densities and electrostatic potentials from file */

  open_now=TRUE;
  for (idim=0;idim<Ndim;idim++) ijk_old_max[idim] = 0;
  for (index=0; index<Nodes_old; index++) {

    if (open_now){
      if( (fp5=fopen(filename,"r")) == NULL){
	printf("Can't open file %s\n", filename);
	exit(1);
      }
       if (Type_poly == CMS || Type_poly==CMS_SCFT || Type_poly==WJDC){
         sprintf(filename2,"%sg",filename);
	 if( (fp6=fopen(filename2,"r")) == NULL){
/*	   printf("Can't open file %s\n", filename2);
	   exit(1);*/
           if ((Iguess_fields=CALC_ALL_FIELDS || Iguess_fields==CALC_RHOBAR_AND_G) && Proc==0) 
              printf("I can't find the file dft_dens.datg... I will construct an initial guess for G equations\n");
           else if ((Iguess_fields=BULK || Iguess_fields==CALC_RHOBAR_ONLY)&& Proc==0) 
              printf("I can't find the file dft_dens.datg... I will construct a simple bulk initial guess for G equations\n");
           Restart_field[G_CHAIN]=FALSE;
	 }
         else  Restart_field[G_CHAIN]=TRUE;
       }
       open_now=FALSE;
                                      /* discard header when ready to read */
       for (i=0;i<header;i++) while ((c=getc(fp5)) != EOF && c !='\n') ; 
    }

    if (Restart==RESTART_1DTOND) ndim_max=1;  /* again for using 1D solution for 2D/3D guess */
    else ndim_max=Ndim;
    /*printf("Proc=%d  ndim_max=%d\n",Proc,ndim_max);*/

                                 /* find number of nodes in each dimension in the file */
    for (idim=0;idim<ndim_max;idim++)         {
      fscanf(fp5,"%lf",&pos_old);
      if (index==0) pos_old0[idim]=pos_old;
      pos_old -= pos_old0[idim];
                                  /* this code transforms x to y etc in the guess */
                                      /*if (idim==0) dim_tmp=1;
                                        else if (idim==1) dim_tmp=0;
                                        else   dim_tmp=idim; */

      dim_tmp=idim;
      ijk_old[dim_tmp] = round_to_int(pos_old/Esize_x[dim_tmp]);
      if ((Type_poly==CMS || Type_poly==CMS_SCFT || Type_poly==WJDC) && Restart_field[G_CHAIN]==TRUE)  fscanf(fp6,"%lf",&tmp); /* ignore positions in densg files. */

      if (ijk_old[dim_tmp] > ijk_old_max[dim_tmp]) ijk_old_max[dim_tmp] = ijk_old[dim_tmp];
    }
                                    /* identify the node and starting unknown number at that node */
    inode=ijk_to_node(ijk_old);
    if (Restart==RESTART_1DTOND) inode=index;
    node_start=inode*Nunk_per_node;

                                   /* loop over unknows assume the order of input is correct */

    for (iunk_file=0;iunk_file<unk_in_file;iunk_file++) {
          eq_type = unk_to_eq_in_file[iunk_file];

          switch (eq_type){
              case CMS_FIELD: 
              case WJDC_FIELD: 
   	         fscanf(fp5,"%lf",&tmp);
                 tmp = exp(-tmp); 
                 break;

              case DENSITY:
                               /* icomp = iunk_file-unk_start_in_file[DENSITY]; */
                               /* if (Lsteady_state) tmp *= Mass[icomp]; */
   	         fscanf(fp5,"%lf",&tmp); 
                 break;

              case HSRHOBAR:
                 if (Restart == RESTART_1DTOND && iunk_file-unk_start_in_file[HSRHOBAR]-Nrho_bar_s>0) tmp=0.0;
                 else fscanf(fp5,"%lf",&tmp); 
                 break;
              case POISSON:
              case CAVWTC:
              case BONDWTC:
   	         fscanf(fp5,"%lf",&tmp); 
                 break;

              case DIFFUSION: 
                               /* icomp = iunk_file-unk_start_in_file[DIFFUSION]; */
                               /* if (Ipot_ff_n != IDEAL_GAS) */  /* Debroglie wavelength contribution */
                               /* tmp -= 3.0*log(Sigma_ff[icomp][icomp]+1.5*log(Mass[icomp]*Temp));*/
   	         fscanf(fp5,"%lf",&tmp); break;
           }
       if (eq_type==DENSITY && convert_to_comp_densities){
            /* this is case where we have read in densities for all the segments and we need to collapse them to component densities */
            icomp=Unk2Comp[iunk_file-unk_start_in_file[eq_type]];
            iunk=Phys2Unk_first[eq_type]+icomp;
            if (Lbinodal && iguess==BINODAL_FLAG) X2_old[iunk+node_start]+=tmp;
            else                                  X_old[iunk+node_start]+=tmp;
       }
       else if (eq_type==DENSITY && convert_to_seg_densities){
            /* this is case where we have read in densities for the components, but we need to expand them to the segment densities -
               note that there is some loss of information here. */
            icomp=iunk_file-unk_start_in_file[eq_type];
            for (jseg=0; jseg<Nseg_tot;jseg++){
               junk=Phys2Unk_first[eq_type]+jseg;
               if (Unk2Comp[jseg]==icomp){
                   if (Lbinodal && iguess==BINODAL_FLAG) X2_old[junk+node_start]+=tmp/Nmer_comp[icomp];
                   else                                  X_old[junk+node_start]+=tmp/Nmer_comp[icomp];
               }
            }
       }
       else{
          iunk = Phys2Unk_first[eq_type]+iunk_file-unk_start_in_file[eq_type];
          if (Lbinodal && iguess==BINODAL_FLAG) X2_old[iunk+node_start]=tmp;
          else                                  X_old[iunk+node_start]=tmp;
       }
    }
 
    if ( (Type_poly==WJDC || Type_poly == CMS || Type_poly==CMS_SCFT) && Restart_field[G_CHAIN]==TRUE){
        for (iunk=Phys2Unk_first[G_CHAIN];iunk<Phys2Unk_last[G_CHAIN];iunk++) {

         /*ipol=Unk_to_Poly[iunk-Phys2Unk_first[G_CHAIN]];  * debugging - CMS ? *
         iseg=Unk_to_Seg[iunk-Phys2Unk_first[G_CHAIN]];
         itype_mer=Phys2Unk_first[CMS_FIELD]+Type_mer[ipol][iseg];*/

         fscanf(fp6,"%lf",&tmp);
         if (Lbinodal && iguess==BINODAL_FLAG) X2_old[iunk+node_start]=tmp;
         else                                  X_old[iunk+node_start]=tmp;

       } /* end of loop over unknowns */
    }  /* end of loading for the CMS G equations. */


             /* read extra variables on the line and ignore them -
                for example the Poisson-Boltzman electrolyte output */
    while ((c=getc(fp5)) != EOF && c !='\n') ;
    if (Restart==RESTART_1DTOND && ijk_old[0]==Nodes_x[0]-1){
       fclose(fp5);
       if ((Type_poly == CMS || Type_poly ==CMS_SCFT || Type_poly==WJDC) && Restart_field[G_CHAIN]==TRUE) fclose(fp6);
       open_now=TRUE;
/* if (Proc==0) printf("closing files to read again!\n");*/
    }
  }
  for (idim=0; idim<Ndim; idim++) {
    Nodes_x_old[idim] = ijk_old_max[idim] + 1;
  }
  if (Restart!=RESTART_1DTOND){
       fclose(fp5);
       if ((Type_poly == CMS || Type_poly ==CMS_SCFT || Type_poly==WJDC) && Restart_field[G_CHAIN]==TRUE) fclose(fp6);
  }
  return;
}
/*******************************************************************/
/*shift_the_profile: do this if the new mesh and the old mesh
                   have identical Esize, but not identical Nnodes_x */
void shift_the_profile(double *x_new,double fac)
{

  int idim,jdim,iunk,inode,inode_old,ijk[3],ijk_tmp[3],Nadd;
  double x_test,unk_old,unk_1,unk_2;

  printf("try to fix up profile in idim=%d in position=%d\n",Plane_new_nodes,Pos_new_nodes);

  idim = Plane_new_nodes;
  Nadd = round_to_int(Del_1[idim]/Esize_x[idim]);
  for (inode=0; inode<Nnodes; inode++){
     node_to_ijk(inode,ijk);
     for (iunk=0; iunk<Nunk_per_node; iunk++){

     switch(Pos_new_nodes){
         case  0:          /*ADDING NODES TO CENTER OF BOX */

           if (ijk[idim] < Nodes_x_old[idim]/2){           /*NODE LBB OF NEW PLANE*/
             inode_old = locate_inode_old(ijk);
             unk_old = X_old[inode_old*Nunk_per_node+iunk];
           }
           else if (ijk[idim] > Nodes_x_old[idim]/2+Nadd){ /*NODE RTF OF NEW PLANE*/
             for (jdim=0; jdim<Ndim; jdim++) ijk_tmp[jdim] = ijk[jdim];
             ijk_tmp[idim] -= Nadd;
             inode_old = locate_inode_old(ijk_tmp);
             unk_old = X_old[inode_old*Nunk_per_node+iunk];
           }
           else {                                    /*NODE IN CENTER OF NEW PLANE*/
             for (jdim=0; jdim<Ndim; jdim++) ijk_tmp[jdim] = ijk[jdim];
             ijk_tmp[idim] = Nodes_x_old[idim]/2;
             inode_old = locate_inode_old(ijk_tmp);
             unk_1 = X_old[inode_old*Nunk_per_node+iunk];

             for (jdim=0; jdim<Ndim; jdim++) ijk_tmp[jdim] = ijk[jdim];
             ijk_tmp[idim] = Nodes_x_old[idim]/2 + 1;
             inode_old = locate_inode_old(ijk_tmp);
             unk_2 = X_old[inode_old*Nunk_per_node+iunk];
      
             unk_old = 0.5*(unk_1+unk_2);
           }
           break;

         case -1:          /*ADDING NODES TO LEFT,BACK,BOTTOM*/
           if (ijk[idim] < Nadd){
             for (jdim=0; jdim<Ndim; jdim++) ijk_tmp[jdim] = ijk[jdim];
             ijk_tmp[idim] = 0;
             inode_old = locate_inode_old(ijk_tmp);
             unk_old = X_old[inode_old*Nunk_per_node+iunk];
           }
           else{
             for (jdim=0; jdim<Ndim; jdim++) ijk_tmp[jdim] = ijk[jdim];
             ijk_tmp[idim] -= Nadd;
             inode_old = locate_inode_old(ijk_tmp);
             unk_old = X_old[inode_old*Nunk_per_node+iunk];
           }
           break;

         case  1:          /*ADDING NODES TO RIGHT,TOP,FRONT*/
           if (ijk[idim] < Nodes_x_old[idim]){
             inode_old = locate_inode_old(ijk);
             unk_old = X_old[inode_old*Nunk_per_node+iunk];
           }
           else {
             for (jdim=0; jdim<Ndim; jdim++) ijk_tmp[jdim] = ijk[jdim];
             ijk_tmp[idim] = Nodes_x_old[idim]-1;
             inode_old = locate_inode_old(ijk_tmp);
             unk_old = X_old[inode_old*Nunk_per_node+iunk];
           }
           break;

         default:
           printf("Check pos new nodes %d [should be -1,0,1]\n",Pos_new_nodes);
           exit(-1);
           break;
 
     }

     x_test=unk_old;
     /* check a few limiting values ... and finally set initial guess*/
/*     if (Unk2Phys[iunk]==DENSITY){
        x_test = AZ_MIN(Rho_max,fac*(unk_old-Rho_b[iunk-Phys2Unk_first[DENSITY]])+ Rho_b[iunk-Phys2Unk_first[DENSITY]]);
     }
     else x_test=unk_old;

     else if (Unk2Phys[iunk]==HSRHOBAR){
        x_test = fac*(unk_old-Rhobar_b[iunk-Phys2Unk_first[HSRHOBAR]])+ Rhobar_b[iunk-Phys2Unk_first[HSRHOBAR]];
        if (iunk == Ncomp && x_test >= 1.0) x_test = Rhobar_b[iunk-Phys2Unk_first[HSRHOBAR]];
     }
     else if (Unk2Phys[iunk]==POISSON){
        x_test = fac*(unk_old-1.0) + 1.0;
     }
     else x_test = unk_old;*/

     x_new[inode*Nunk_per_node+iunk] = x_test;

     }
  }
  return;
}
/*********************************************************************/
/* locate_inode_old: find a given node in the old solution array */
int locate_inode_old(int *ijk)
{
   int inode_old=0;
 
   if (Ndim == 1)      inode_old = ijk[0];
   else if (Ndim == 2) inode_old = ijk[0] + ijk[1] * Nodes_x_old[0];
   else if (Ndim == 3) inode_old = ijk[0] + ijk[1] * Nodes_x_old[0] 
                       + ijk[2] * Nodes_x_old[0]*Nodes_x_old[1];
   return(inode_old);
}
/**********************************************************************/
/* communicate profile: broadcast the x_new profile to all processors,
                        and let each of them pick out the needed entries.*/
void communicate_profile(double *x_new, double** xInBox)
{
   int inode,iunk,inode_box;   
  
    MPI_Bcast (x_new, Nnodes*Nunk_per_node,MPI_DOUBLE,0,MPI_COMM_WORLD);

    for (inode_box=0; inode_box<Nnodes_box; inode_box++){
       inode = B2G_node[inode_box];
       for (iunk=0; iunk<Nunk_per_node; iunk++){
           xInBox[iunk][inode_box] = x_new[inode*Nunk_per_node+iunk];
       }
    }
    return;
}
/*********************************************************************/
/*check_zero_densities: here just remove zero densities where 
         not appropriate, and make sure x=0 where needed */
void check_zero_densities(double **xInBox)
{

  int loc_inode,icomp,inode_box,iunk,iloop,nloop;

  if (Lseg_densities) nloop=Nseg_tot;
  else nloop=Ncomp;

  for (inode_box=0; inode_box<Nnodes_box; inode_box++){
      for (iloop=0; iloop<nloop; iloop++){
         icomp=Unk2Comp[iloop];
	 iunk = Phys2Unk_first[DENSITY]+iloop;
         if (Zero_density_TF[inode_box][icomp])
                 xInBox[iunk][inode_box] = 0.0;
         else{
           if (Lseg_densities)
              if (xInBox[iunk][inode_box] < Rho_seg_b[iunk]*exp(-VEXT_MAX)) {
                  xInBox[iunk][inode_box] = Rho_seg_b[iunk]*exp(-VEXT_MAX); /*DENSITY_MIN*/
              }
           else
              if (xInBox[iunk][inode_box] < Rho_b[icomp]*exp(-VEXT_MAX)) {
                  xInBox[iunk][inode_box] = Rho_b[icomp]*exp(-VEXT_MAX); /*DENSITY_MIN*/
              }
         }
      }
  }
  return;
}
/**********************************************************************/
/*chop_profile: do the profile chop here. */
void chop_profile(double **xInBox, int iguess)
{
  int loc_inode,icomp,iwall,check,iunk,inode_box;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++)
    inode_box=L2B_node[loc_inode];
    for (icomp=0; icomp<Ncomp; icomp++) {
        check = 0;
        for (iwall=0; iwall<Nwall; iwall++)
              if (X_wall[loc_inode][iwall] > Xstart_step[0]) check++;
        if (check == Nwall){
	      iunk = Phys2Unk_first[DENSITY]+icomp;
              if (iguess==CHOP_RHO_L) xInBox[iunk][inode_box] = Rho_coex[1];
              else if (iguess==CHOP_RHO_V) xInBox[iunk][inode_box] = Rho_coex[0];
              else if (iguess==CHOP_RHO) xInBox[iunk][inode_box] = Rho_b[icomp];
              else if (iguess==CHOP_RHO_STEP) xInBox[iunk][inode_box] = Rho_step[icomp][0];
        }
    }
  return;
}
/****************************************************************************/
/*setup_exp_density_with_profile: in this routine set up a density
                     profile as rho(x)*exp(-Vext/kT)*/
void setup_exp_density_with_profile(double **xInBox)
{

  int loc_inode,i,inode_box,iunk,icomp;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     for (i=0;i<Ncomp;i++){
        iunk = icomp+Phys2Unk_first[DENSITY];
        if (Vext[loc_inode][icomp]>0.0) xInBox[iunk][inode_box] *= exp(-Vext[loc_inode][i]);

/*        if (Type_poly==CMS || Type_poly==CMS_SCFT)
        xInBox[i+Phys2Unk_first[CMS_FIELD]][inode_box] = exp(-log(xInBox[i+Phys2Unk_first[CMS_FIELD]][inode_box])+Vext[loc_inode][i]);*/
     }
  }

  return;
}
/****************************************************************************/

