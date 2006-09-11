/*******NOTE - THIS BIT OF CODE SHOULD BE PUT UP ONE LEVEL IN THE PKG DIRECTORY, BUT IT SHOULD BE GENERALIZED 
        WITH FNC POINTERS AND GENERIC STRINGS .....if I can figure out how to do it..... ************/


/* fill_HSRHOBAR_eqn: top level routine for filling the nonlocal denstity (rhobar)  equations
void fill_HSRHOBAR_eqn(inode,iunk,resid_only_flag)
{
  int ieq;

   for(ieq=0;ieq<NEQ_TYPE;ieq++){
      if (Block_flag[ieq] != ZERO_BLOCK_FLAG){
      if (Phys2Nunk[ieq] !=0) {
         switch(ieq){
                                                  /* nonzero matrix block */
            case DENSITY:     Block_flag[ieq]=fill_HSRHOBAR_DENSITY(inode,iunk,resid_only_flag); break;

                                                  /* identity matrix block */
            case HSRHOBAR:    Block_flag[ieq]= fill_HSRHOBAR_HSRHOBAR(inode,iunk,resid_only_flag); break;

                                                  /* zero matrix blocks - see below */
            case POISSON:     Block_flag[ieq]= fill_HSRHOBAR_POISSON(inode,iunk,resid_only_flag); break;
            case DIFFUSION:   Block_flag[ieq]= fill_HSRHOBAR_DIFFUSION(inode,iunk,resid_only_flag); break;
            case CAVITY_WTC:  Block_flag[ieq]= fill_HSRHOBAR_CAVITY_WTC(inode,iunk,resid_only_flag); break;
            case BOND_WTC:    Block_flag[ieq]= fill_HSRHOBAR_BOND_WTC(inode,iunk,resid_only_flag); break;
            case USR_VAR1:    Block_flag[ieq]= fill_HSRHOBAR_USR_VAR1(inode,iunk,resid_only_flag); break;
            case USR_VAR2:    Block_flag[ieq]= fill_HSRHOBAR_USR_VAR2(inode,iunk,resid_only_flag); break;
            case USR_VAR3:    Block_flag[ieq]= fill_HSRHOBAR_USR_VAR3(inode,iunk,resid_only_flag); break;
         }
       }
     }
   }
   return;
}
