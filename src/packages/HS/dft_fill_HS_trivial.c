/*****************************************************************************/
int fill_HSRHOBAR_HSRHOBAR(int inode,int iunk,int resid_only_flag)
{
  resid =-x[iunk][inode_box];
  resid_sum+=resid;
  mat_val=-1.0;
  dft_linprobmgr_insertrhsvalue(LinProbMgr_manager,iunk,loc_inode,-resid);
  dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,mat_val);

  return (IDENTITY_BLOCK_FLAG);
}
/*****************************************************************************/
int fill_HSRHOBAR_POISSON(int inode,int iunk,int resid_only_flag)
{
   return (ZERO_BLOCK_FLAG);
}
/*****************************************************************************/
int fill_HSRHOBAR_DIFFUSION(int inode,int iunk,int resid_only_flag)
{
   return (ZERO_BLOCK_FLAG);
}
/*****************************************************************************/
int fill_HSRHOBAR_CAVITY_WTC(int inode,int iunk,int resid_only_flag)
{
   return (ZERO_BLOCK_FLAG);
}
/*****************************************************************************/
int fill_HSRHOBAR_USR_VAR1(int inode,int iunk,int resid_only_flag)
{
   return (ZERO_BLOCK_FLAG);
}
/*****************************************************************************/
int fill_HSRHOBAR_USR_VAR2(int inode,int iunk,int resid_only_flag)
{
   return (ZERO_BLOCK_FLAG);
}
/*****************************************************************************/
int fill_HSRHOBAR_USR_VAR3(int inode,int iunk,int resid_only_flag)
{
   return (ZERO_BLOCK_FLAG);
}

