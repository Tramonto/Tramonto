#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

/****************************************************************************/
double integrand_elec_PB_freen(int iunk,int inode_box, double **x)
{
     double integrand,rho_i,psi_i;
     int icomp,psiunk;

     if (Type_poly==WTC) icomp = Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     else icomp = iunk-Phys2Unk_first[DENSITY];
     
     rho_i = x[iunk][inode_box];
     psiunk = Phys2Unk_first[POISSON];
     psi_i = x[psiunk][inode_box];

     
     if (rho_i > 0.) integrand = 0.5*rho_i*Charge_f[icomp]*psi_i;
     return(integrand);
}
/****************************************************************************/
double integrand_elec_MSAcorr_freen(int iunk,int inode_box, double **x)
{
     double integrand,rho_i;

     rho_i = x[iunk][inode_box];

     int_stencil(x,inode_box,iunk,THETA_CHARGE);
     if (rho_i > 0.0) integrand = -0.5*rho_i*Temporary_sum;
     return(integrand);
}
/****************************************************************************/
double integrand_elec_MSAcorr_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_bulk;
     int icomp,iseg;

     if (Type_poly==WTC){
         iseg = iunk-Phys2Unk_first[DENSITY];
         icomp = Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
         if(Lsteady_state) rho_bulk = Rho_seg_RTF[iseg];
         else              rho_bulk = Rho_seg_b[iseg];
     }
     else{
           icomp = iunk-Phys2Unk_first[DENSITY];
           if(Lsteady_state) rho_bulk = Rho_b_RTF[icomp];
           else              rho_bulk = Rho_b[icomp];
     }

     integrand = -0.5*rho_bulk*Deltac_b[icomp];
     return(integrand);
}
/****************************************************************************/

