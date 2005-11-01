#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

/****************************************************************************/
double integrand_vext_freen(int iunk,int inode_box, double **x)
{     double integrand,rho_i;
     int icomp,iseg,i;

     i = iunk-Phys2Unk_first[DENSITY];
     if (Type_poly==WTC) icomp = Unk2Comp[i];
     else                icomp = i;

     rho_i = x[iunk][inode_box];

     if (rho_i > DENSITY_MIN) integrand = rho_i*Vext[B2L_node[inode_box]][icomp];
     
     return(integrand);
}
/****************************************************************************/
double integrand_vext_elec_freen(int iunk,int inode_box, double **x)
{     double integrand,rho_i;
     int icomp,iseg,i;

     i = iunk-Phys2Unk_first[DENSITY];
     if (Type_poly==WTC) icomp = Unk2Comp[i];
     else                icomp = i;

     rho_i = x[iunk][inode_box];

     if (rho_i > DENSITY_MIN) integrand = 0.5*rho_i*Vext_coul[B2L_node[inode_box]][icomp];
     
     return(integrand);
}
