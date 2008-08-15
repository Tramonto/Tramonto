/*
 *  dft_guess_SCF.c
 *  tramonto_cassia
 *
 *  Created by alfrisc on 8/8/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "dft_guess_SCF.h"

/*********************************************************/
/*setup_polymer_field: in this routine set up the initial guess for the CMS-SCF field variable */
/* called if not calculating all fields; just a guess for the fields */
void setup_polymer_SCF_field(double **xInBox, int iguess)
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
		}
	}
	return;
}
/*********************************************************/
/*calc_init_SCFfield: calculate the SCF field from knowledge of rho and lambda if available */
void calc_init_SCFfield(double **xInBox)
{
	int loc_inode,icomp,jcomp,jrho,iunk,unk_L,inode_box;
	double field,int_bulk;
	
	printf("in calc_init_SCFfield\n");
	
	for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
		inode_box=L2B_node[loc_inode];
		
		unk_L = Phys2Unk_first[SCF_CONSTR];
		
		for (icomp=0; icomp<Ncomp; icomp++){
			iunk = Phys2Unk_first[SCF_FIELD]+icomp;
			if (!Zero_density_TF[inode_box][icomp]){
				field= Vext[loc_inode][icomp];
				if(Restart_field[SCF_CONSTR])
					field += xInBox[unk_L][inode_box];
				for (jcomp=0; jcomp<Ncomp; jcomp++){
					if(jcomp != icomp) {
						jrho = Phys2Unk_first[DENSITY]+jcomp;
						field -= Eps_ff[icomp][jcomp]*xInBox[jrho][inode_box]/Rho_t;
					}
				}
			}
			else field=VEXT_MAX;
			xInBox[iunk][inode_box]=exp(-field);		
		}
	}
	return;
}

/*********************************************************/
/*calc_init_lambda: set up the initial guess for the lambda constraint variable */
void calc_init_lambda(double **xInBox)
{
	int loc_inode,inode_box,iunk;
	
	printf("in calc_init_lambda\n");
	
	iunk = Phys2Unk_first[SCF_CONSTR];
	
	for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
		inode_box=L2B_node[loc_inode];
		xInBox[iunk][inode_box]=0.0;
	}
	return;
}


