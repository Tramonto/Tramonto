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
 *  FILE: dft_guess_SCF.c
 *
 *  Created by alfrisc on 8/8/08.
 *  This file contains routines to set up an initial guess for SCF theory.
 *
 */

#include "dft_guess_SCF.h"

/*********************************************************/
/*setup_polymer_field: in this routine set up the initial guess for the CMS-SCF field variable */
/* called if not calculating all fields; just a guess for the fields */
void setup_polymer_SCF_field(double **xInBox, double **xOwned, int iguess)
{
	int loc_inode,itype_mer,irho, iunk,inode_box;
	double field;
	printf("in setup_polymer_SCF_field\n");
	
	for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
		inode_box=L2B_node[loc_inode];
		
		for (itype_mer=0; itype_mer<Ncomp; itype_mer++){
			irho = Phys2Unk_first[DENSITY]+itype_mer;
			iunk = Phys2Unk_first[SCF_FIELD]+itype_mer;
			if (xInBox[irho][inode_box]<1.e-6) field=VEXT_MAX-1.;
			else field=-log(xInBox[irho][inode_box]/Rho_b[itype_mer]);
			xInBox[iunk][inode_box]=exp(-field);
                        xOwned[iunk][loc_inode]=xInBox[iunk][inode_box];
		}
	}
	return;
}
/*********************************************************/
/*calc_init_SCFfield: calculate the SCF field from knowledge of rho and lambda if available */
void calc_init_SCFfield(double **xInBox, double **xOwned)
{
	int loc_inode,icomp,jcomp,jrho,iunk,unk_L,inode_box;
	double field,int_bulk;
	
	field = 0.0;

	
	for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
		inode_box=L2B_node[loc_inode];
		
		for (icomp=0; icomp<Ncomp; icomp++){
			iunk = Phys2Unk_first[SCF_FIELD]+icomp;
			if (!Zero_density_TF[inode_box][icomp]){
				field= Vext[loc_inode][icomp];
				if(Phys2Nunk[SCF_CONSTR] && Restart_field[SCF_CONSTR]) {
					unk_L = Phys2Unk_first[SCF_CONSTR];
					field += xInBox[unk_L][inode_box];
				}
				for (jcomp=0; jcomp<Ncomp; jcomp++){
					jrho = Phys2Unk_first[DENSITY]+jcomp;
					field += Eps_ff[icomp][jcomp]*xInBox[jrho][inode_box]/Rho_t;
				}
			}
			else field=VEXT_MAX;
			xInBox[iunk][inode_box]=exp(-field);		
                        xOwned[iunk][loc_inode]=xInBox[iunk][inode_box];
		}
	}
	return;
}

/*********************************************************/
/*calc_init_lambda: set up the initial guess for the lambda constraint variable */
void calc_init_lambda(double **xInBox,double **xOwned)
{
	int loc_inode,inode_box,iunk;
	
	iunk = Phys2Unk_first[SCF_CONSTR];
	
	for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
		inode_box=L2B_node[loc_inode];
		xInBox[iunk][inode_box]=0.0;
                xOwned[iunk][loc_inode]=xInBox[iunk][inode_box];
	}
	return;
}

/*********************************************************/
/*calc_init_polymer_G_SCF: in this routine sets up the initial guess for the chain variable
in the SCF functional 
this routine uses the machinery in load_Chain_Geqns_SCF and particularly load_polymer_recursion to do all
the integrals to calculate the G's from the initial conditions and values of rho and the fields */
void calc_init_polymer_G_SCF(double **xInBox,double **xOwned)
{
	int loc_inode,itype_mer,irho, iunk,i,Nloop,inode_box,field;
	int ibond,jbond,index,iseg,jseg,pol_num,bond_num,test,ijk_box[3];
	double resid_G;
	int array_val[NMER_MAX*NBOND_MAX],array_fill,count_fill;
	double (*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **);
	double (*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **);
	
	fp_ResidG=&CMS_Resid_GCHAIN;
	fp_ResidG_Bulk=&CMS_Resid_Bulk_GCHAIN;
	
	/* need to be careful to generate the G's in the order dictated
		by the chain architecture.  Use same strategy as in dft_thermo_wjdc */
	
	field=SCF_FIELD;
	
	for (ibond=0;ibond<Nbonds;ibond++) array_val[ibond]=FALSE;
	array_fill=FALSE;
	count_fill=0;
	
	while (array_fill==FALSE){
		for (ibond=0;ibond<Nbonds;ibond++){
			pol_num=Unk_to_Poly[ibond];
			iseg=Unk_to_Seg[ibond];
			bond_num=Unk_to_Bond[ibond];
			
			if (array_val[ibond]==FALSE){
				test=TRUE;  /* assume we will compute a bulk G */
				jseg=Bonds[pol_num][iseg][bond_num];
				if (jseg != -1 && jseg != -2){   /* may need to skip this G if we don't have all information yet  -
					always compute G for end segments flagged with -1 value */
					for (jbond=0;jbond<Nbond[pol_num][jseg];jbond++){
						if (Bonds[pol_num][jseg][jbond] != iseg){ /* check all jbonds to see if we have necessary info */
							index=Poly_to_Unk[pol_num][jseg][jbond]+Geqn_start[pol_num]-Geqn_start[0];
							if (array_val[index]==FALSE) test=FALSE;
						}
					}
				}
				if (test==TRUE){     /* compute a bulk G */
                                       (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xOwned, xInBox);  /* make sure all densities are available for calculations */

					for (loc_inode=0;loc_inode<Nnodes_per_proc;loc_inode++){
						inode_box=L2B_node[loc_inode];
						node_box_to_ijk_box(inode_box,ijk_box);
						iunk=Phys2Unk_first[G_CHAIN]+ibond;
						resid_G=load_Chain_Geqns_SCF(field,0,0,
												 NULL,fp_ResidG,fp_ResidG_Bulk,
												 iunk,loc_inode,inode_box,
												 ijk_box,0,xInBox, INIT_GUESS_FLAG);
						
						xInBox[iunk][inode_box]=resid_G;
                                                xOwned[iunk][loc_inode]=xInBox[iunk][inode_box];
					}
					count_fill++;
					array_val[ibond]=TRUE;
					
				} /*end of test=TRUE condition */
			}
		}
		if (count_fill==Nbonds) array_fill=TRUE;
	}
return;
}
/*********************************************************/


