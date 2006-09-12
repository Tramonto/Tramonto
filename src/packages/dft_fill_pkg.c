/****************************************************************************/

/* top level routine - dft_fill_main.c */

in if (Unk2Phys[]==XXX) loop of dft_fill_main.c we will need to change logic to read (for example)
if (Unk2Phys[iunk]==HSRHOBAR)
   fill_package(eq_hsrhobar,inode,iunk,resid_only_flag,x,resid);
etc...

we will also need to add a routine to make the function pointer assignments early on in the code - it only needs to be done
once - the heart of the routine should be as follows where for each set of equations, function pointer assignments are made.
This is a bit ugly, but it will avoid ugliness in the fill routines later. 

/* function pointer definitions for EULER-LAGRANGE EQUATIONS */

/* function pointer definitions for HSRHOBAR EQUATIONS */

/* function pointer definitions for POISSONS EQUATION */
eq_poisson.(*density_fill)=&fill_poisson_density;
eq_poisson.(*hsrhobar_fill)=&fill_poisson_hsrhobar;
eq_poisson.(*poisson_fill)=&fill_poisson_poisson;
eq_poisson.(*diffusion_fill)=&fill_poisson_diffusion;
eq_poisson.(*cavwtc_fill)=&fill_poisson_cavwtc;
eq_poisson.(*bondwtc_fill)=&fill_poisson_bondwtc;
eq_poisson.(*cmsfield_fill)=&fill_poisson_cmsfield;
eq_poisson.(*cmsgeqn_fill)=&fill_poisson_cmsgeqn;
eq_poisson.(*usrvar1_fill)=&fill_poisson_usrvar1;
eq_poisson.(*usrvar2_fill)=&fill_poisson_usrvar2;
eq_poisson.(*usrvar3_fill)=&fill_poisson_usrvar3;

/* function pointer definitions for DIFFUSION EQUATIONS */

/* function pointer definitions for WTC CAVITY  EQUATIONS */

/* function pointer definitions for WTC BOND EQUATIONS */

/* function pointer definitions for CMS FIELD EQUATIONS */

/* function pointer definitions for CMS G EQUATIONS */

/* function pointer definitions for USERFN1 EQUATIONS */

/* function pointer definitions for USERFN2 EQUATIONS */

/* function pointer definitions for USERFN3 EQUATIONS */


need to 
1. assign function addresses to the various fill_phys structures
2. send the correct structure to the fill_package routine below.

/****************************************************************************/
/* fill_package: routine that takes a fill_phys struct, puts it into a generic name, eq_type, and then
   calls the necessary associated functions to fill a given row of the matrix */
double fill_package(fill_phys eq_type,int loc_inode,int iunk,int icomp,int inode_box, int izone, int *ijk_box,
                     int resid_only_flag,double **x,double *resid)
{
  int iphys,flag;

   for(iphys=0;iphys<NEQ_TYPE;iphys++){
      if (Phys2Nunk[iphys] !=0) {
      if (Block_flag[Unk2Phys[iunk]][iphys] != ZERO_BLOCK_FLAG){
         switch(iphys){
            case DENSITY:     
                 flag=eq_type->(*density_fill)(loc_inode,inode_box,iunk,icomp,inode_box,izone,ijk_box,resid_only_flag,x,resid);
                 break;
            case HSRHOBAR:    
                 flag=eq_type->(*hsrhobar_fill)(loc_inode,inode_box,iunk,icomp,inode_box,izone,ijk_box,resid_only_flag,x,resid);
                 break;
            case POISSON:     
                 flag=eq_type->(*poisson_fill)(loc_inode,inode_box,iunk,icomp,inode_box,izone,ijk_box,resid_only_flag,x,resid);
                 break;
            case DIFFUSION:   
                 flag=eq_type->(*diffusion_fill)(loc_inode,inode_box,iunk,icomp,inode_box,izone,ijk_box,resid_only_flag,x,resid);
                 break;
            case CAVITY_WTC:  
                 flag=eq_type->(*cavwtc_fill)(loc_inode,inode_box,iunk,icomp,inode_box,izone,ijk_box,resid_only_flag,x,resid);
                 break;
            case BOND_WTC:    
                 flag=eq_type->(*bondwtc_fill)(loc_inode,inode_box,iunk,icomp,inode_box,izone,ijk_box,resid_only_flag,x,resid);
                 break;
            case USR_VAR1:    
                 flag=eq_type->(*usrvar1_fill)(loc_inode,inode_box,iunk,icomp,inode_box,izone,ijk_box,resid_only_flag,x,resid);
                 break;
            case USR_VAR2:    
                 flag=eq_type->(*usrvar2_fill)(loc_inode,inode_box,iunk,icomp,inode_box,izone,ijk_box,resid_only_flag,x,resid);
                 break;
            case USR_VAR3:    
                 flag=eq_type->(*usrvar3_fill)(loc_inode,inode_box,iunk,icomp,inode_box,izone,ijk_box,resid_only_flag,x,resid);
                 break;
         }
         Block_flag[Unk2Phys[iunk]][iphys]=flag;
       }
     }
   }
   return (resid);
}
/*****************************************************************************/
