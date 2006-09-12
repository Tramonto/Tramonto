/****************************************************************************/
/* setup_package: routine that sets up the function pointers based on the desired physics in the calculation of interest */
double setup_package()
{
  int iphys,flag;

   for(iphys=0;iphys<NEQ_TYPE;iphys++){
      if (Phys2Nunk[iphys] !=0) {
         switch(iphys){
            case DENSITY:     
                 eq_el.(*density_fill)=&fill_el_density;
                 eq_el.(*hsrhobar_fill)=&fill_el_hsrhobar;
                 eq_el.(*poisson_fill)=&fill_el_poisson;
                 eq_el.(*diffusion_fill)=&fill_el_diffusion;
                 eq_el.(*cavwtc_fill)=&fill_el_cavwtc;
                 eq_el.(*bondwtc_fill)=&fill_el_bondwtc;
                 eq_el.(*cmsfield_fill)=&fill_el_cmsfield;
                 eq_el.(*cmsgeqn_fill)=&fill_el_cmsgeqn;
                 eq_el.(*usrvar1_fill)=&fill_el_usrvar1;
                 eq_el.(*usrvar2_fill)=&fill_el_usrvar2;
                 eq_el.(*usrvar3_fill)=&fill_el_usrvar3;
                 break;

            case HSRHOBAR:    
                 eq_hsrhobar.(*density_fill)=&fill_hsrhobar_density;
                 eq_hsrhobar.(*hsrhobar_fill)=&fill_hsrhobar_hsrhobar;
                 eq_hsrhobar.(*poisson_fill)=&fill_hsrhobar_poisson;
                 eq_hsrhobar.(*diffusion_fill)=&fill_hsrhobar_diffusion;
                 eq_hsrhobar.(*cavwtc_fill)=&fill_hsrhobar_cavwtc;
                 eq_hsrhobar.(*bondwtc_fill)=&fill_hsrhobar_bondwtc;
                 eq_hsrhobar.(*cmsfield_fill)=&fill_hsrhobar_cmsfield;
                 eq_hsrhobar.(*cmsgeqn_fill)=&fill_hsrhobar_cmsgeqn;
                 eq_hsrhobar.(*usrvar1_fill)=&fill_hsrhobar_usrvar1;
                 eq_hsrhobar.(*usrvar2_fill)=&fill_hsrhobar_usrvar2;
                 eq_hsrhobar.(*usrvar3_fill)=&fill_hsrhobar_usrvar3;
                 break;

            case POISSON:     
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
                 break;

            case DIFFUSION:   
                 eq_diffusion.(*density_fill)=&fill_diffusion_density;
                 eq_diffusion.(*hsrhobar_fill)=&fill_diffusion_hsrhobar;
                 eq_diffusion.(*poisson_fill)=&fill_diffusion_poisson;
                 eq_diffusion.(*diffusion_fill)=&fill_diffusion_diffusion;
                 eq_diffusion.(*cavwtc_fill)=&fill_diffusion_cavwtc;
                 eq_diffusion.(*bondwtc_fill)=&fill_diffusion_bondwtc;
                 eq_diffusion.(*cmsfield_fill)=&fill_diffusion_cmsfield;
                 eq_diffusion.(*cmsgeqn_fill)=&fill_diffusion_cmsgeqn;
                 eq_diffusion.(*usrvar1_fill)=&fill_diffusion_usrvar1;
                 eq_diffusion.(*usrvar2_fill)=&fill_diffusion_usrvar2;
                 eq_diffusion.(*usrvar3_fill)=&fill_diffusion_usrvar3;
                 break;

            case CAVWTC:  
                 eq_cavwtc.(*density_fill)=&fill_cavwtc_density;
                 eq_cavwtc.(*hsrhobar_fill)=&fill_cavwtc_hsrhobar;
                 eq_cavwtc.(*poisson_fill)=&fill_cavwtc_poisson;
                 eq_cavwtc.(*diffusion_fill)=&fill_cavwtc_diffusion;
                 eq_cavwtc.(*cavwtc_fill)=&fill_cavwtc_cavwtc;
                 eq_cavwtc.(*bondwtc_fill)=&fill_cavwtc_bondwtc;
                 eq_cavwtc.(*cmsfield_fill)=&fill_cavwtc_cmsfield;
                 eq_cavwtc.(*cmsgeqn_fill)=&fill_cavwtc_cmsgeqn;
                 eq_cavwtc.(*usrvar1_fill)=&fill_cavwtc_usrvar1;
                 eq_cavwtc.(*usrvar2_fill)=&fill_cavwtc_usrvar2;
                 eq_cavwtc.(*usrvar3_fill)=&fill_cavwtc_usrvar3;
                 break;

            case BONDWTC:    
                 eq_bondwtc.(*density_fill)=&fill_bondwtc_density;
                 eq_bondwtc.(*hsrhobar_fill)=&fill_bondwtc_hsrhobar;
                 eq_bondwtc.(*poisson_fill)=&fill_bondwtc_poisson;
                 eq_bondwtc.(*diffusion_fill)=&fill_bondwtc_diffusion;
                 eq_bondwtc.(*cavwtc_fill)=&fill_bondwtc_cavwtc;
                 eq_bondwtc.(*bondwtc_fill)=&fill_bondwtc_bondwtc;
                 eq_bondwtc.(*cmsfield_fill)=&fill_bondwtc_cmsfield;
                 eq_bondwtc.(*cmsgeqn_fill)=&fill_bondwtc_cmsgeqn;
                 eq_bondwtc.(*usrvar1_fill)=&fill_bondwtc_usrvar1;
                 eq_bondwtc.(*usrvar2_fill)=&fill_bondwtc_usrvar2;
                 eq_bondwtc.(*usrvar3_fill)=&fill_bondwtc_usrvar3;
                 break;

            case USRFNC1:    
                 eq_usrfnc1.(*density_fill)=&fill_usrfnc1_density;
                 eq_usrfnc1.(*hsrhobar_fill)=&fill_usrfnc1_hsrhobar;
                 eq_usrfnc1.(*poisson_fill)=&fill_usrfnc1_poisson;
                 eq_usrfnc1.(*diffusion_fill)=&fill_usrfnc1_diffusion;
                 eq_usrfnc1.(*cavwtc_fill)=&fill_usrfnc1_cavwtc;
                 eq_usrfnc1.(*bondwtc_fill)=&fill_usrfnc1_bondwtc;
                 eq_usrfnc1.(*cmsfield_fill)=&fill_usrfnc1_cmsfield;
                 eq_usrfnc1.(*cmsgeqn_fill)=&fill_usrfnc1_cmsgeqn;
                 eq_usrfnc1.(*usrvar1_fill)=&fill_usrfnc1_usrvar1;
                 eq_usrfnc1.(*usrvar2_fill)=&fill_usrfnc1_usrvar2;
                 eq_usrfnc1.(*usrvar3_fill)=&fill_usrfnc1_usrvar3;
                 break;
            case USRFNC2:    
                 eq_usrfnc2.(*density_fill)=&fill_usrfnc2_density;
                 eq_usrfnc2.(*hsrhobar_fill)=&fill_usrfnc2_hsrhobar;
                 eq_usrfnc2.(*poisson_fill)=&fill_usrfnc2_poisson;
                 eq_usrfnc2.(*diffusion_fill)=&fill_usrfnc2_diffusion;
                 eq_usrfnc2.(*cavwtc_fill)=&fill_usrfnc2_cavwtc;
                 eq_usrfnc2.(*bondwtc_fill)=&fill_usrfnc2_bondwtc;
                 eq_usrfnc2.(*cmsfield_fill)=&fill_usrfnc2_cmsfield;
                 eq_usrfnc2.(*cmsgeqn_fill)=&fill_usrfnc2_cmsgeqn;
                 eq_usrfnc2.(*usrvar1_fill)=&fill_usrfnc2_usrvar1;
                 eq_usrfnc2.(*usrvar2_fill)=&fill_usrfnc2_usrvar2;
                 eq_usrfnc2.(*usrvar3_fill)=&fill_usrfnc2_usrvar3;
                 break;
            case USRFNC3:    
                 eq_usrfnc3.(*density_fill)=&fill_usrfnc3_density;
                 eq_usrfnc3.(*hsrhobar_fill)=&fill_usrfnc3_hsrhobar;
                 eq_usrfnc3.(*poisson_fill)=&fill_usrfnc3_poisson;
                 eq_usrfnc3.(*diffusion_fill)=&fill_usrfnc3_diffusion;
                 eq_usrfnc3.(*cavwtc_fill)=&fill_usrfnc3_cavwtc;
                 eq_usrfnc3.(*bondwtc_fill)=&fill_usrfnc3_bondwtc;
                 eq_usrfnc3.(*cmsfield_fill)=&fill_usrfnc3_cmsfield;
                 eq_usrfnc3.(*cmsgeqn_fill)=&fill_usrfnc3_cmsgeqn;
                 eq_usrfnc3.(*usrvar1_fill)=&fill_usrfnc3_usrvar1;
                 eq_usrfnc3.(*usrvar2_fill)=&fill_usrfnc3_usrvar2;
                 eq_usrfnc3.(*usrvar3_fill)=&fill_usrfnc3_usrvar3;
                 break;
         }
     }
   }
   return (resid);
}
/*****************************************************************************/
