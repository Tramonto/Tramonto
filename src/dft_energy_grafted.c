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

#include "dft_energy_grafted.h"

/* calculate chain contributions to free energy of mixture of grafted and non-grafted polymers --
 essentially convert grand free energy to a semi-grand free energy by Legendre transform for 
 grafted polymers */

double WJDCgraft_freen(double **x)
{
    int i,pol_comp,icomp,iunk,inode;
    double omega[NCOMP_MAX],fhelm[NCOMP_MAX],fid[NCOMP_MAX];
    double om[NCOMP_MAX][NCOMP_MAX],fh[NCOMP_MAX][NCOMP_MAX];
    double omega_wjdc,fidtot,ftot,freen_tot,rho_i;
    
    omega_wjdc = 0.0;
    fidtot = 0.0;
    ftot = 0.0;

    /* loop over all polymers */
    for(pol_comp=0; pol_comp<Npol_comp; pol_comp++) {
        omega[pol_comp] = 0.0;
        fhelm[pol_comp] = 0.0;

        /* loop over components in polymer */
        for(icomp=0;icomp<Ncomp;icomp++){
            iunk = Phys2Unk_first[DENSITY]+icomp;
            if(pol_comp==Icomp_to_polID[icomp]) {
            if (Lseg_densities)
                om[pol_comp][icomp]=integrateInSpace(&integrand_WJDC_freen,iunk,Nel_hit2,x,NULL);
            else
                om[pol_comp][icomp]=integrateInSpace(&integrand_WJDCcomp_freen,iunk,Nel_hit2,x,NULL);
             omega[pol_comp] += om[pol_comp][icomp];
            }
        }
        if (Proc==0 && Iwrite_screen == SCREEN_VERBOSE) printf("ipol=%d,omega=%f\n",pol_comp,omega[pol_comp]);
        if(Grafted[pol_comp]) {
            /* calculate extra term in semi-grand free energy for grafted chains */
            /* loop over components in polymer */
            for(icomp=0;icomp<Ncomp;icomp++){
                iunk = Phys2Unk_first[DENSITY]+icomp;
                if(pol_comp==Icomp_to_polID[icomp]) {
                fh[pol_comp][icomp] = graft_freen(pol_comp,iunk,icomp,x);
                fhelm[pol_comp] += fh[pol_comp][icomp];
                }
            }
            if (Proc==0 && Iwrite_screen == SCREEN_VERBOSE) printf("ipol=%d,fh=%f\n",pol_comp,fhelm[pol_comp]);
        }
        omega_wjdc += omega[pol_comp];
        ftot += fhelm[pol_comp];
    }
    if (Proc==0 && Iwrite_screen == SCREEN_VERBOSE) print_to_screen(omega_wjdc,"WJDC term");
    if (Proc==0 && Iwrite_screen == SCREEN_VERBOSE) print_to_screen(ftot,"WJDC Helmholtz term, grafted");
    freen_tot = omega_wjdc + ftot;
    return(freen_tot);
}

    
/****************************************************************************/
/* calculate the term integral of mu_i*rho where mu_i = chain segment chemical potential */
double graft_freen(int npol, int iunk, int icomp, double **x)
{
    double sum_i,sum,rho_lam,lambda;
    int loc_inode,inode_box;
    
    sum_i = 0.0;
    sum = 0.0;
    if(Grafted[npol]==GRAFT_DENSITY)
        lambda = log(Rho_g[npol]/Gsum_graft[npol]);
    else
        lambda = log(Rho_g[npol]/(Gsum_graft[npol]*Total_area_graft[npol]));
    if (Proc==0 && Iwrite_screen == SCREEN_VERBOSE) print_to_screen(lambda,"grafted chain chem. potential");
    /* integrate over mu*rho here */
    for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
        inode_box = L2B_node[loc_inode];
        
        /*divide by chain length to get segment chemical potential */
        rho_lam = lambda*x[iunk][inode_box]/Nmer[npol];
     /*   if(Grafted_TypeID[npol]==icomp){
            sum_i += rho_lam;
        }
        else    */
            sum_i += rho_lam*Nel_hit2[icomp][inode_box]*Vol_el/((double)Nnodes_per_el_V);
    }
    sum=gsum_double(sum_i);
    
    return(sum);
}

