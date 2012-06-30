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
 *  FILE: dft_fill_SCF.c
 *
 *  Created by alfrisc on 6/9/08.
 *
 * This file contains the matrix fill routine associated with the SCF theory for polymers.
 */

#include "dft_fill_SCF.h"


/****************************************************************************/
double load_SCF_density(int iunk, int loc_inode, int inode_box, double **x,int resid_only_flag) 
{
	int itype_mer,unk_B;
	double resid_R;
		
	resid_R=0.0;
	
	itype_mer = iunk-Phys2Unk_first[DENSITY];
		
	if (Zero_density_TF[inode_box][itype_mer] || Vext[loc_inode][itype_mer] == VEXT_MAX){
		resid_R=fill_zero_value(iunk,loc_inode,inode_box,x,resid_only_flag);
	}
	else{
		unk_B=Phys2Unk_first[SCF_FIELD]+itype_mer; /* Boltzmann factor for this i */
		resid_R+=resid_and_Jac_ChainDensity (G_CHAIN,x,iunk,unk_B,loc_inode,inode_box,
							 resid_only_flag, &prefactor_rho_scft);
	}
	return(resid_R);
}

/****************************************************************************/
double prefactor_rho_scft(int itype_mer,int inode_box,double **x)
{
	int npol;
	double fac;

	npol = 0;
	while (Nmer_t[npol][itype_mer]==0) npol++;
	fac = Rho_b[itype_mer]/Nmer_t[npol][itype_mer];
		
	return (fac);
}


/****************************************************************************/
double load_SCF_field(int iunk, int loc_inode, int inode_box, int *ijk_box, int izone, double **x,int resid_only_flag) 
{
	double resid_B,mat_val;
	int itype_mer,junk,unk_L,jcomp;
		
	itype_mer = iunk - Phys2Unk_first[SCF_FIELD];
	
    if (Zero_density_TF[inode_box][itype_mer] || Vext[loc_inode][itype_mer] == VEXT_MAX){
		resid_B=fill_zero_value(iunk,loc_inode,inode_box,x,resid_only_flag);
    }
    else{
		resid_B=0.0;
		for(jcomp=0; jcomp<Ncomp; jcomp++) {
			junk = Phys2Unk_first[DENSITY] + jcomp;
			resid_B += Eps_ff[itype_mer][jcomp]*x[junk][inode_box]/Rho_t;
		}
		if(Phys2Nunk[SCF_CONSTR]) {
			unk_L=Phys2Unk_first[SCF_CONSTR];  /* lambda field */
			resid_B += x[unk_L][inode_box];
		}
		resid_B += Vext[loc_inode][itype_mer]+log(x[iunk][inode_box]);
		if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
			dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid_B);
		if(!resid_only_flag){
			mat_val = 1.0/x[iunk][inode_box];  /* delta R / delta field */
			dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);
			for(jcomp=0; jcomp<Ncomp; jcomp++) {   /* delta R / delta rho */
				junk = Phys2Unk_first[DENSITY] + jcomp;
				mat_val = Eps_ff[itype_mer][jcomp]/Rho_t;
				dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,junk,inode_box,mat_val);
			}
			if(Phys2Nunk[SCF_CONSTR]) {
				mat_val = 1.0;	/* delta R / delta lambda */
				dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_L,inode_box,mat_val);
			}
		}		
    }
    return(resid_B);
}

/****************************************************************************/
/* incompressibility constraint for SCF */
double load_lambda_field(int iunk, int loc_inode, int inode_box, int *ijk_box, int izone, double **x,int resid_only_flag) 
{
	int icomp,unk_rho;
	double resid_L, mat_val;
	
	resid_L = 0.0;

	for(icomp=0; icomp<Ncomp; icomp++){
		/* decide how to treat surfaces! */
		if(Zero_density_TF[inode_box][icomp] || Vext[loc_inode][icomp] == VEXT_MAX){
			resid_L=fill_zero_value(iunk,loc_inode,inode_box,x,resid_only_flag);
		}
		else{
			unk_rho=Phys2Unk_first[DENSITY]+icomp;
			resid_L += x[unk_rho][inode_box];
			if(!resid_only_flag){
				mat_val = 1.0;
				dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,unk_rho,inode_box,mat_val);
			}
		}	
	}
	
	resid_L -= Rho_t;
	
	if(!Zero_density_TF[inode_box][icomp] && Vext[loc_inode][icomp] != VEXT_MAX) {
		if (resid_only_flag != INIT_GUESS_FLAG && resid_only_flag != CALC_RESID_ONLY) 
			dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid_L);
	}
	
	return(resid_L);
}

/****************************************************************************/
double load_SCF_Geqns(int iunk, int loc_inode, int inode_box, int *ijk_box, int izone, double **x,int resid_only_flag)
{
    int Njacobian_types;
    int Njacobian_sums;
	int field;
    double resid_G;
    void (*funcArray_Jac[3])(int,int,int,int,int,int,int,int,int *,double,double **);
    double (*fp_ResidG)(int,int,int,int,int,int,int,int *,double,double **);
    double (*fp_ResidG_Bulk)(int,int,int,int,int,int,int,int *,double,double **);
	
    Njacobian_types=2;   /* derivatives with respect to field and G are needed */
    Njacobian_sums=1;    /* a sum(integral) must be computed for the field derivatives */
	
    funcArray_Jac[0]=&CMS_Jacobian_GCHAIN_derivG;
    funcArray_Jac[1]=&CMS_Jacobian_GCHAIN_derivFIELD;
    fp_ResidG=&CMS_Resid_GCHAIN;
    fp_ResidG_Bulk=&CMS_Resid_Bulk_GCHAIN;
	
	field=SCF_FIELD;
	
    resid_G=load_Chain_Geqns_SCF(field,Njacobian_types,Njacobian_sums,
                             funcArray_Jac,fp_ResidG,fp_ResidG_Bulk,
                             iunk, loc_inode,inode_box, 
                             ijk_box,izone,x, resid_only_flag);
    return(resid_G);
}

/*******************************************************************************
collect_G_old: This gathers the G that we need on Proc 0.        */ 

/*void collect_G_old(double **x, int unk_G, double *G_old)
{
	int i,iunk,loc_inode, loc_i,idim,nunk_per_proc;
	int *index=NULL;
	double *unk_global, *unk_loc;
	
	Nodes_old = Nnodes;
	for (idim=0; idim<Ndim; idim++) Nodes_x_old[idim] = Nodes_x[idim];*/
	
	/* allocate temporary arrays */
	
/*	unk_loc = (double *) array_alloc (1, Nnodes_per_proc, sizeof(double));
	
	for (loc_inode=0; loc_inode < Nnodes_per_proc; loc_inode++ ) */
/*		unk_loc[loc_inode] = x[unk_G][L2B_node[loc_inode]];  * always use nodal ordering here */
					
/*	if (Proc == 0) {
		unk_global = (double *) array_alloc (1, Nnodes, sizeof(double));
		index = (int *) array_alloc (1, Nnodes, sizeof(int));
	}
	else {
		unk_global=NULL;
		index=NULL;
	}*/
			
	/* collect the node numbers from all the processors */
			
/*	MPI_Gatherv(L2G_node,Nnodes_per_proc,MPI_INT,
				index,Comm_node_proc,Comm_offset_node,
				MPI_INT,0,MPI_COMM_WORLD);*/
	
	/* collect the unknowns from all the processors */
	
/*	MPI_Gatherv(unk_loc,Nnodes_per_proc,MPI_DOUBLE,
				unk_global,Comm_unk_proc,Comm_offset_unk,
				MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	safe_free((void *) &unk_loc);
	
	if (Proc == 0){
		for (i=0; i<Nnodes; i++)*/
/*			G_old[index[i]] = unk_global[i];  * figure out correct indexing */
/*		safe_free((void *) &unk_global);
		safe_free((void *) &index);
	}
	
	return;
}*/

/****************************************************************************/
/*double load_SCF_density(int iunk, int loc_inode, int inode_box, double **x,int resid_only_flag) 
{
	int itype_mer, Ns;
	double resid;

	itype_mer = iunk-Phys2Unk_first[DENSITY];
	
	resid = 0.0;
	return(resid);
}*/

/****************************************************************************/
/* allocate memory for the q propagators */
/*void setup_SCF_q() {
	double delta_s;
	int Ns;
	
	delta_s = 1.0;
	Ns = Nmer[0]/delta_s;*/

	/* create propagator functions; q[icomp][s][r] */
/*	q = (double ***) array_alloc (3, Ncomp, Ns+1, Nnodes_box, sizeof(double));

	return;
}*/


/****************************************************************************/
/* function to solve the diffusion equation for the chain propagators */
/*void solve_q_diff(double **x) {
	
	int s, loc_inode, inode_box;
	int unk_w, Ns;
	double q_tmp, delta_s; */
	/* static variables keep their value for every time the function is called*/
/*	static double *wt_lp_1el, *wt_s_1el;
	static int   **elem_permute;	
	
	Ns = Nmer[0];
	delta_s = 1.0;*/
	
	/* for now,  set field equal to a function */
/*	unk_w = Phys2Unk_first[SCF_FIELD];
	for  (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {*/
		
		/* convert local node to box */
/*		inode_box = L2B_node[loc_inode];
		x[unk_w][inode_box] = 1.0;
	}*/
	
	/* load weights for Laplacian */
/*	set_fem_1el_weights(&wt_lp_1el, &wt_s_1el, &elem_permute);*/
	
	/* loop over spatial nodes */
	/* will have to deal with ghost nodes for q somehow!! */
/*	for  (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++) {
		inode_box = L2B_node[loc_inode];*/
/*		q[0][0][inode_box] = 1.0;		* initial conditions for q */
/*		for(s=0; s<Ns; s++){			* loop along chain */
/*			q_tmp = exp(-delta_s*x[unk_w][inode_box]/2.0)*q[0][s][inode_box];
			q[0][s+1][inode_box] = q_tmp;
		}
	}
	
	
	
	return;
}*/

