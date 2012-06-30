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
 *  FILE: dft_fill.c
 *
 *  This file contains the fill of the residual equations and Jacobian
 *  matrix.
 */

#include "dft_fill_main.h"

/****************************************************************************/
double fill_resid_and_matrix (double **x, struct RB_Struct *dphi_drb, int iter, int resid_only_flag,int unk_flag)
{
 /*
  * Local variable declarations
  */

  char   *yo = "fill_resid_and_matrix";
  int     loc_inode, inode_box,ijk_box[3],iunk,iunk_start,iunk_end;
  int     mesh_coarsen_flag_i;
  /* int switch_constmatrix;*/
  int	npol,iseg,unk_G,idim,inode;
  double nodepos[3];
  double *resid_unk,resid_sum=0.0,resid_term;
  double sum_i, vol;
  double **array_test;
  double Ads[NCOMP_MAX],pos_xyz[3];
  FILE *fparray;
  char filename[FILENAME_LENGTH];
  char *fileArray;
  char tmp_str_array[FILENAME_LENGTH];

  int i,j,ipol,iseg_test,icomp;

  if (Proc == 0 && !resid_only_flag && Iwrite_screen == SCREEN_VERBOSE) printf("\n\t%s: Doing fill of residual and matrix\n",yo);
  resid_unk = (double *) array_alloc (1, Nunk_per_node, sizeof(double));


  if (Iwrite_files==FILES_DEBUG_MATRIX){
      Array_test = (double **) array_alloc (2, Nunk_per_node*Nnodes,Nunk_per_node*Nnodes, sizeof(double));
      for (i=0;i<Nunk_per_node*Nnodes;i++){
          for (j=0;j<Nunk_per_node*Nnodes;j++){ Array_test[i][j]=0.0; }}
  }

  /* for debugging print out profiles on each iteration */
  if (Iwrite_files==FILES_DEBUG) print_profile_box(x, "dens_iter.dat");

  if (unk_flag == NODAL_FLAG){
      iunk_start = 0;
      iunk_end = Nunk_per_node;
  } 
  else{
      iunk_start = unk_flag;
      iunk_end = unk_flag+1;
  }
  
  /* ALF: calculate all the single chain partition functions */
  if(Type_poly == CMS_SCFT) {
	  
	  vol = 1.0;
	  for(idim=0; idim<Ndim; idim++)
		vol *= Size_x[idim];
	  
	  /* loop over chains */
	  for(npol=0; npol<Npol_comp; npol++){
		  iseg = Nmer[npol]-1;
		  unk_G = Geqn_start[npol] + Poly_to_Unk[npol][iseg][0];		/* find G_{N-1}^{N-2} */
		  sum_i=0.0, Gsum[npol] = 0.0;
		  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
			  inode_box = L2B_node[loc_inode];
			  /* need to figure out correct index for Nel_hit2: 1st index is over icomp */
			  /* integrate by extended trapezoidal rule */
			  sum_i += x[unk_G][inode_box]*Nel_hit2[0][inode_box]*Vol_el/((double)Nnodes_per_el_V);
		  }
		 Gsum[npol] = gsum_double(sum_i)/vol; 
		  if(Gsum[npol] < 1.e-6) if (Iwrite_screen==SCREEN_VERBOSE) printf("Proc=%d error: Gsum small = %f\n",Proc,Gsum[npol]);
	  }
  } 
 
	/* ALF: for grafted chains */
/*	if(Type_poly==CMS || Type_poly==WJDC3) calc_Gsum(x);*/
	if(Type_poly==WJDC3 && Grafted_Logical==TRUE){ calc_Gsum_new(x); }
		
		/* Load residuals and matrix */
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {

    /* convert local node to global */

    inode_box = L2B_node[loc_inode];
    node_box_to_ijk_box(inode_box, ijk_box);
    if ( ((Mesh_coarsening != FALSE) && (Nwall_type >0)) || L1D_bc) mesh_coarsen_flag_i = Mesh_coarsen_flag[inode_box];
    else mesh_coarsen_flag_i = 0;
    for (iunk=iunk_start; iunk<iunk_end; iunk++) {
      resid_term=0.0;

      resid_unk[iunk]=0.0;

      if (mesh_coarsen_flag_i == FLAG_1DBC) load_coarse_node_1dim(loc_inode,inode_box,ijk_box,iunk,x,resid_only_flag);

      else if (mesh_coarsen_flag_i < 0 && 
               mesh_coarsen_flag_i != FLAG_BULK && 
               mesh_coarsen_flag_i != FLAG_PBELEC) load_coarse_node_Ndim(loc_inode,inode_box,iunk,x,resid_only_flag);

      else{
          /*switch_constmatrix=FALSE;
          if (iter>1 && resid_only_flag==FALSE && Constant_row_flag[Unk2Phys[iunk]]==TRUE) {
             resid_only_flag=CALC_RESID_ONLY; switch_constmatrix=TRUE;
          }*/
          resid_term=load_standard_node(loc_inode,inode_box,ijk_box,iunk,x,dphi_drb,
                               resid_unk,mesh_coarsen_flag_i,resid_only_flag);
          resid_sum+=resid_term;
/*          if (switch_constmatrix) resid_only_flag=FALSE;*/
      }

     /* print for debugging purposes call this print routine */ 
       if (Iwrite_screen==SCREEN_DEBUG_RESID)  print_residuals(loc_inode,iunk,resid_unk);

    } /* end of loop over # of unknowns per node */
  } /* end of loop over local nodes */

  if (Iwrite_files==FILES_DEBUG_MATRIX){
     sprintf(filename,"matrix_proc%0d.dat",Proc);
     fileArray=filename;
     strcpy(tmp_str_array,Outpath_array);
     fparray=fopen(strcat(tmp_str_array,fileArray),"w");

     for (i=0;i<Nunk_per_node*Nnodes;i++){
        for (j=0;j<Nunk_per_node*Nnodes;j++){
            if (fabs(Array_test[i][j])>1.e-12) fprintf(fparray,"%d  %d  %d  %g\n",Proc,Nunk_per_node*Nnodes-j,i,Array_test[j][i]);
     }}
     safe_free((void *) &Array_test);
     fclose(fparray);
  }

  if(Type_poly==WJDC3 && Grafted_Logical==TRUE){
           safe_free((void *) &Index_SurfNodes_Gsum);
           safe_free((void *) &Index_UnkGQ_Gsum);
    /*       safe_free((void *) &Gsum_graft);
           safe_free((void *) &Gsum_graft_noVolume);*/
           safe_free((void *) &Nodes_Surf_Gsum);
           safe_free((void *) &Total_area_graft);
           safe_free((void *) &GsumPrefac_XiDerivs);
           safe_free((void *) &GsumPrefac_GDerivs);
  }

  safe_free((void *) &resid_unk);
/*  if (Type_poly==CMS_SCFT) safe_free((void *) &Gsum);*/
  return(resid_sum);
}

/*************************************************************************************************/
void calc_Gsum_new(double **x)
{
  double *integrand,*integrand2;
  double prefac,prefac2,gproduct,gproduct_deriv;
  int  surf_norm,idim,iwall_type_graft,iel_w,jbond;
  int icomp,ilist,ipol,ibond,unk_GQ,loc_inode,iseg_graft,isegALL_graft,iwall,unk_B,Nbond_graft,inode,inode_box;

  integrand = (double *) array_alloc (1, Nwall, sizeof(double));
  integrand2 = (double *) array_alloc (1, Nwall, sizeof(double));

  if(Type_poly==WJDC3 && Grafted_Logical ==TRUE) {
    Nodes_Surf_Gsum = (int *) array_alloc (1, Nwall, sizeof(int));
    Index_SurfNodes_Gsum = (int **) array_alloc(2, Nwall, Nnodes_per_proc, sizeof(int));
    if (Gsum_graft==NULL){
    Gsum_graft = (double *) array_alloc(1, Npol_comp, sizeof(double));
    Gsum_graft_noVolume = (double *) array_alloc(1, Npol_comp, sizeof(double));
    }
    Total_area_graft = (double *) array_alloc(1, Npol_comp, sizeof(double));

    Index_UnkGQ_Gsum = (int ***) array_alloc(3, Nwall, Nnodes_per_proc, NBOND_MAX, sizeof(int));
    Index_UnkB_Gsum = (int **) array_alloc(2, Nwall, Nnodes_per_proc, sizeof(int));
    GsumPrefac_XiDerivs = (double **) array_alloc (2, Nwall, Nnodes_per_proc, sizeof(double));
    GsumPrefac_GDerivs = (double ***) array_alloc (3, Nwall, Nnodes_per_proc, NBOND_MAX, sizeof(double));

    for(ipol=0; ipol<Npol_comp; ipol++){
      Gsum_graft[ipol]=1.0;
      if(Grafted[ipol]) {

        iseg_graft=Grafted_SegID[ipol];
        isegALL_graft=Grafted_SegIDAll[ipol];
        iwall_type_graft=Graft_wall[ipol];
        icomp=Grafted_TypeID[ipol];
        Nbond_graft=Nbonds_SegAll[isegALL_graft];
        unk_B=Phys2Unk_first[WJDC_FIELD]+icomp;

        if (Nlists_HW == 1 || Nlists_HW == 2) ilist = 0;
        else                                  ilist = icomp;

        for (iwall=0;iwall<Nwall;iwall++){
            integrand[iwall]=0.0;
            integrand2[iwall]=0.0;
            Nodes_Surf_Gsum[iwall]=0;
        }

/*        for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) 
            if (Nodes_2_boundary_wall[ilist][L2B_node[loc_inode]] !=-1)*/

        for (inode=0;inode<NodesS_global[ilist];inode++){
          inode_box=S2B_node[ilist][inode];

/*        for (iel_w=0; iel_w<Nelems_S[ilist][loc_inode]; iel_w++)
            surf_norm = Surf_normal[ilist][loc_inode][iel_w];*/

        for (iel_w=0; iel_w<NelemsS_global[ilist][inode]; iel_w++){
            surf_norm = Surf_normal_global[ilist][inode][iel_w];
 
            idim = abs(surf_norm) - 1;

/*            iwall=Surf_elem_to_wall[ilist][loc_inode][iel_w];*/
            iwall=Surf_elem_to_wall_global[ilist][inode][iel_w];

            if (WallType[iwall]==iwall_type_graft){
    /*           prefac = Area_surf_el[idim]*Esize_x[idim]/Nelems_S[ilist][loc_inode];
               prefac2 = Area_surf_el[idim]/Nelems_S[ilist][loc_inode];*/

               prefac = Area_surf_el[idim]*Esize_x[idim]/NelemsS_global[ilist][inode];
               prefac2 = Area_surf_el[idim]/NelemsS_global[ilist][inode];

               Index_UnkB_Gsum[iwall][Nodes_Surf_Gsum[iwall]]=unk_B;
               GsumPrefac_XiDerivs[iwall][Nodes_Surf_Gsum[iwall]]=
                    ((double)(Nbond_graft-2))*prefac/(POW_DOUBLE_INT(x[unk_B][inode_box],Nbond_graft-1));
 
               prefac/=POW_DOUBLE_INT(x[unk_B][inode_box],Nbond_graft-2);
               prefac2/=POW_DOUBLE_INT(x[unk_B][inode_box],Nbond_graft-2);
 
               gproduct=1.0; 
               for(ibond=0; ibond<Nbonds_SegAll[isegALL_graft]; ibond++) {
                    if(Bonds[ipol][iseg_graft][ibond] != -1) {

                       unk_GQ = Geqn_start[ipol] + Poly_to_Unk[ipol][iseg_graft][ibond];
                       gproduct*=x[unk_GQ][inode_box];

                       gproduct_deriv=1.0;
                       for (jbond=0; jbond<Nbonds_SegAll[isegALL_graft]; jbond++){
                          if(Bonds[ipol][iseg_graft][jbond] != -1 && (Bonds[ipol][iseg_graft][ibond] != Bonds[ipol][iseg_graft][jbond])) {
                          unk_GQ = Geqn_start[ipol] + Poly_to_Unk[ipol][iseg_graft][jbond];
                          gproduct_deriv*=x[unk_GQ][inode_box];
                          }
                       }
                       Index_UnkGQ_Gsum[iwall][Nodes_Surf_Gsum[iwall]][ibond]=unk_GQ;
                       GsumPrefac_GDerivs[iwall][Nodes_Surf_Gsum[iwall]][ibond]=gproduct_deriv*prefac;
                    }
               }
               integrand[iwall]+=prefac*gproduct;
               integrand2[iwall]+=prefac2*gproduct;
               GsumPrefac_XiDerivs[iwall][Nodes_Surf_Gsum[iwall]]*=gproduct;
               Index_SurfNodes_Gsum[iwall][Nodes_Surf_Gsum[iwall]]=S2B_node[ilist][inode];
               Nodes_Surf_Gsum[iwall]+=1;
            }
         }
         }

         Gsum_graft[ipol]=0.0;
         Total_area_graft[ipol]=0.0;
         Gsum_graft_noVolume[ipol]=0.0;
         for (iwall=0;iwall<Nwall;iwall++){
            if (WallType[iwall]==iwall_type_graft){
               Gsum_graft[ipol]+=(integrand[iwall]);
/*               Gsum_graft_noVolume[ipol]+=2.0*integrand2[iwall];*/
               Total_area_graft[ipol]+=S_area_tot[ilist][iwall];
            }
         }
         Gsum_graft[ipol]/=Total_area_graft[ipol];
/*         Gsum_graft_noVolume[ipol]/=Total_area_graft[ipol];*/
        }
     }
  }

  safe_free((void *) &integrand);
  safe_free((void *) &integrand2);

  return;
}
/*************************************************************************************************/

void calc_Gsum(double **x) 
{
	int npol,ibond,iseg,itype_mer,unk_G;
	int graft_seg,graft_bond,gbond,izone, unk_xi2,unk_xi3;
	int iwall,pwall,pwall_type,inode_w;
	int loc_inode,inode,inode_box,ijk_box[3];
	double gsum,nodepos[3],nodeposc[3];
	double y,ysqrt,xi_2,xi_3;	
	
	for(npol=0; npol<Npol_comp; npol++){
		gsum = 0.0;
		if(Grafted[npol]) {
			/* find type of bonded site */
			for(iseg=0; iseg<Nmer[npol]; iseg++) {
				for(ibond=0; ibond<Nbonds_SegAll[iseg]; ibond++) {
					if(Bonds[npol][iseg][ibond] == -2) {	/* grafted site */
						itype_mer = Type_mer[npol][iseg];
						graft_seg = iseg;
						graft_bond = ibond;
					}
				}
			}			
			/* find correct G; assumes 2 bonds only for grafting site */
			for (ibond=0; ibond<Nbonds_SegAll[graft_seg]; ibond++)
				if(ibond != graft_bond) gbond = ibond;
			unk_G = Geqn_start[npol] + Poly_to_Unk[npol][graft_seg][gbond];
			 
			 /* find wall position and do integral with stencil */
			for(loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {
				inode = L2B_node[loc_inode];
				inode_box = L2B_node[loc_inode];
				iwall = Nodes_2_boundary_wall[0][inode_box];
				if(iwall != -1) {
					node_to_position(inode,nodepos); 
					
					/* correct wall? */
					pwall = Graft_wall[npol];
					pwall_type = WallType[pwall];
					nodeposc[0] = WallPos[0][pwall] + Poly_graft_dist[pwall_type];
					nodeposc[1] = nodeposc[2] = 0.0;
					if(nodepos[0]==nodeposc[0] + Sigma_ff[itype_mer][itype_mer]/2.0 || 
					   nodepos[0]==nodeposc[0] - Sigma_ff[itype_mer][itype_mer]/2.0) {
						inode_w = position_to_node(nodeposc);
						inode_box = node_to_node_box(inode_w);
					node_box_to_ijk_box(inode_box, ijk_box);
					if ( ((Mesh_coarsening != FALSE) && (Nwall_type >0)) || L1D_bc) izone = Mesh_coarsen_flag[inode_box];
					else izone = 0;
					if(Type_poly==CMS) {
						gsum = grafted_int(DELTA_FN_BOND,itype_mer,ijk_box,izone,unk_G,x);
					}
					else if(Type_poly == WJDC3) {
						/* find y at wall position */
						unk_xi2=Phys2Unk_first[CAVWTC]; unk_xi3=Phys2Unk_first[CAVWTC]+1;
						xi_2=x[unk_xi2][inode_box]; xi_3=x[unk_xi3][inode_box];
						y=y_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[itype_mer][itype_mer],xi_2,xi_3);
						ysqrt=sqrt(y); 
						gsum = ysqrt*grafted_int(DELTA_FN_BOND,itype_mer,ijk_box,izone,unk_G,x); 
					}
					else {
						if(Proc==0 && Iwrite_screen != SCREEN_NONE) printf("grafting not implemented for Type_poly=%d\n", Type_poly);
						exit(-1);
					}
					}					
				}
			} /* end loop over loc_inode */
			
			/* communicate value of Gsum */
			Gsum[npol] = gsum_double(gsum);
		}
		else
			Gsum[npol] = 1.0;
	}
	
	return;
}


/*************************************************************************************************/
/* grafted_int:  this function calculates an integral for overall density normalization with grafted chains */
double grafted_int(int sten_type,int itype_mer, int *ijk_box, int izone, int unk_G,double **x)
{
	int   **sten_offset, *offset, isten;
	int unk_xi2,unk_xi3;
	double *sten_weight,  weight;
	struct Stencil_Struct *sten;
	double y,ysqrt,xi_2,xi_3;
	
	int jnode_box,jcomp;
	int reflect_flag[NDIM_MAX];
	double sum = 0.0;
	
	sten = &(Stencil[sten_type][izone][itype_mer+Ncomp*itype_mer]);
	sten_offset = sten->Offset;
	sten_weight = sten->Weight;
	
	for (isten = 0; isten < sten->Length; isten++) {
		offset = sten_offset[isten];
		weight = sten_weight[isten];
		
		/* Find the Stencil point */
		jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);
		jcomp=itype_mer;
		if (jnode_box >= 0 && !Zero_density_TF[jnode_box][jcomp]) {			
			if(Type_poly == CMS)
				sum += weight*x[unk_G][jnode_box];
			else if (Type_poly == WJDC3) {
				unk_xi2=Phys2Unk_first[CAVWTC]; unk_xi3=Phys2Unk_first[CAVWTC]+1;
				xi_2=x[unk_xi2][jnode_box]; xi_3=x[unk_xi3][jnode_box];
				y=y_cav(Sigma_ff[itype_mer][itype_mer],Sigma_ff[itype_mer][itype_mer],xi_2,xi_3);
				ysqrt=sqrt(y);
				sum += weight*ysqrt*x[unk_G][jnode_box];
			}
		}
	}
			
	return(sum);
}


/*************************************************************************************************/
/* load_standard_node:  this routine controls the invocation of different physics routines that load the 
                        matrix problem of interest.  */
double load_standard_node(int loc_inode,int inode_box, int *ijk_box, int iunk, double **x,
                        struct  RB_Struct *dphi_drb, double *resid_unk,
                        int mesh_coarsen_flag_i,int resid_only_flag)
{
   int izone,icomp;
   double l2norm_term;
                                           
   /* IZONE: izone is the zone number that controls the quadrature scheme to be used. */

   izone = mesh_coarsen_flag_i;

   switch(Unk2Phys[iunk]){
       case DENSITY: 
             if (L_HSperturbation && Type_poly != WJDC && Type_poly != WJDC2 && Type_poly != WJDC3){
                resid_unk[iunk]=load_euler_lagrange(iunk,loc_inode,inode_box,ijk_box,izone,x,
                                                    dphi_drb,mesh_coarsen_flag_i,resid_only_flag);
             }
             else if(Type_poly==CMS  || Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3) {                              
                if (Type_poly==CMS) resid_unk[iunk]=load_CMS_density(iunk,loc_inode,inode_box,x,resid_only_flag);
                else                resid_unk[iunk]=load_WJDC_density(iunk,loc_inode,inode_box,x,resid_only_flag);
             }
	     else if(Type_poly == CMS_SCFT) {
		 resid_unk[iunk]=load_SCF_density(iunk,loc_inode,inode_box,x,resid_only_flag); 
	     }
	     else if(Type_poly == SCFT) {
	          if (Iwrite_screen != SCREEN_NONE) printf("SCFT not yet implemented\n");
		  exit(-1);
	     }
             break;
       case HSRHOBAR: 
          if (iunk == Phys2Unk_first[HSRHOBAR]){
             resid_unk[iunk]=load_rho_bar_s(THETA_FN_R,x,iunk,loc_inode,inode_box,izone,ijk_box, resid_only_flag);
          }
          else if (iunk < Phys2Unk_first[HSRHOBAR]+Nrho_bar_s){
             resid_unk[iunk]=load_rho_bar_s(DELTA_FN_R,x,iunk,loc_inode,inode_box,izone,ijk_box, resid_only_flag);
          }
          else if (iunk >= Phys2Unk_first[HSRHOBAR]+Nrho_bar_s){
              resid_unk[iunk]=load_rho_bar_v(x,iunk,loc_inode,inode_box,izone,ijk_box, resid_only_flag);
          }
        break;

       case MF_EQ:
           icomp=iunk-Phys2Unk_first[MF_EQ];
           if (Type_poly != CMS && Type_poly != CMS_SCFT){
                resid_unk[iunk]=load_mean_field(THETA_PAIRPOT_RCUT,iunk,loc_inode,
                                          icomp,izone,ijk_box, x, resid_only_flag);
           }
           else{
                resid_unk[iunk]=load_mean_field(THETA_CR_DATA,iunk,loc_inode,
                                          icomp,izone,ijk_box,x,resid_only_flag);
           }
           break;

       case POISSON: 
		resid_unk[iunk]=load_poisson_control(iunk,loc_inode,inode_box,ijk_box,x,resid_only_flag); break;

       case DIFFUSION:
           if (Linear_transport==TRUE)
               resid_unk[iunk]=load_linear_transport_eqn(iunk,loc_inode,inode_box,ijk_box,x,resid_only_flag);
           else
               resid_unk[iunk]=load_nonlinear_transport_eqn(iunk,loc_inode,inode_box,ijk_box,x,resid_only_flag);
           break;

       case CAVWTC: 
          resid_unk[iunk]=load_cavity_wtc(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
          break;

       case BONDWTC:
          resid_unk[iunk]=load_bond_wtc(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
          break;

       case CMS_FIELD: 
           resid_unk[iunk]=load_CMS_field(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
           break;

       case WJDC_FIELD: 
          resid_unk[iunk]=load_WJDC_field(iunk,loc_inode,inode_box,ijk_box,izone,x,dphi_drb,mesh_coarsen_flag_i,resid_only_flag);
          break;

       case G_CHAIN:
          if (Type_poly==CMS){
             resid_unk[iunk]=load_CMS_Geqns(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
          }
		  else if (Type_poly==CMS_SCFT)
			  resid_unk[iunk]=load_SCF_Geqns(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
          else if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){
             resid_unk[iunk]=load_WJDC_Geqns(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
          }
          break;
       case SCF_FIELD:
	  if(Type_poly==CMS_SCFT) resid_unk[iunk]=load_SCF_field(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
	  else{
	      if (Proc==0 && Iwrite_screen != SCREEN_NONE) printf("SCFT not yet implemented\n");
	      exit(-1);
	  }
	  break;
       case SCF_CONSTR:
	   resid_unk[iunk]=load_lambda_field(iunk,loc_inode,inode_box,ijk_box,izone,x,resid_only_flag);
	   break; 

   }  /* end of physics switch */

   l2norm_term = resid_unk[iunk];  
   l2norm_term *= l2norm_term;

   return(l2norm_term);    /*note that the l2norm term is used for convergence tests - picard iters */
}
/*************************************************************************************************/
/* print_residuals: tool to assist in debugging residual fills.  This routine prints the
   residuals at a given node for each unknown in the problem. Note that these tools are 
   most often used for debugging physics on a single processor, and may need modification
   for multi-processor debugging. */
   void print_residuals(int loc_inode,int iunk,double *resid_unk)
{
               /* note: translate local coordinates to global coordinates with L2G_node[loc_inode]*/
               /* note: translate local coordinates to box coordinates with L2B_node[loc_inode]*/
               /* note: if you want to print to a file you need to modify print statements below */

   char filename[FILENAME_LENGTH]="./resid.out";
               /* also note: to separate parts of the physics constributions (e.g. different parts 
                  of the euler-lagrange equation you will need to modify the return parameters from the
                  various load_ functions. Note that this
                  kind of analysis may require multiple runs and so output to a file is recommended. */

    /* PRINT STATEMENTS FOR PHYSICS DEBUGGING .... CHECK RESIDUALS INDEPENDENTLY  */
/*    if (fabs(resid_unk[iunk])>1.e-3){*/
    switch(Unk2Phys[iunk]){
       case DENSITY:  printf("Proc=%d: loc_inode=%d of %d (Global val=%d) iunk_rho=%d ", Proc,loc_inode,Nnodes_per_proc,L2G_node[loc_inode],iunk); break;
       case HSRHOBAR: printf("Proc=%d: loc_inode=%d iunk_rbar=%d ", Proc,loc_inode,iunk); break;
       case POISSON:  printf("Proc=%d: loc_inode=%d iunk_poisson=%d ", Proc,loc_inode,iunk); break;
       case DIFFUSION: printf("Proc=%d: loc_inode=%d  iunk_diffusion=%d ",Proc,loc_inode,iunk); break;
       case CAVWTC: printf("Proc=%d: loc_inode=%d  iunk_cavity=%d ",Proc,loc_inode,iunk); break;
       case BONDWTC: printf("Proc=%d: loc_inode=%d  iunk_bondwtc=%d ",Proc,loc_inode,iunk); break;
       case CMS_FIELD: printf("Proc=%d: loc_inode=%d  iunk_cmsfield=%d ",Proc,loc_inode,iunk); break;
	   case SCF_FIELD: printf("Proc=%d: loc_inode=%d  iunk_cmsfield=%d ",Proc,loc_inode,iunk); break;
	   case SCF_CONSTR: printf("Proc=%d: loc_inode=%d  iunk_cmsfield=%d ",Proc,loc_inode,iunk); break;
       case WJDC_FIELD: printf("Proc=%d: loc_inode=%d  iunk_wjdc_field=%d ",Proc,loc_inode,iunk); break;
       case G_CHAIN: printf("Proc=%d: loc_inode=%d  iunk_Gchain=%d ",Proc,loc_inode,iunk); break;
       case MF_EQ: printf("Proc=%d: loc_inode=%d  iunk_MFeq=%d ",Proc,loc_inode,iunk); break;
    }
    printf(" resid=%11.8f \n",resid_unk[iunk]); 
/*    }*/

    return;
}
/*************************************************************************************************/
