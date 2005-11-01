#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

/****************************************************************************/
double integrand_att_freen(int iunk,int inode_box, double **x)
{
     double integrand=0.0,rho_i;
     rho_i = x[iunk][inode_box];

     int_stencil(x,inode_box,iunk,U_ATTRACT);
     if (rho_i > 0.) integrand = 0.5*rho_i*Temporary_sum;
     return(integrand);
}
/****************************************************************************/
double integrand_att_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_bulk;
     int icomp,i;

     i= iunk-Phys2Unk_first[DENSITY];
     if (Type_poly==WTC){
           icomp = Unk2Comp[i];
           if (Lsteady_state) rho_bulk = Rho_seg_RTF[i];
           else rho_bulk = rho_bulk = Rho_seg_b[i];
     }
     else{
           icomp = i;
           if (Lsteady_state) rho_bulk = Rho_b_RTF[i];
           else               rho_bulk = Rho_b[i];
     }

     integrand = 0.5*rho_bulk*Betamu_att[icomp];
     return(integrand);
}
/****************************************************************************/
