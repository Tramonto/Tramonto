using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_OutputParams_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Output_List)
{
   string str_screenout,str_fileout,str_printtooutput,str_adsout,str_energyout,str_LperArea,str_reflect,str_pmf;

  /*******************************/
  /* Define a few useful strings */
  /*******************************/
  str_screenout="Select how much output to the screen you would like to generate.\n Basic output includes thermodynamics and post processing results, but minimal output from solver.  \nVerbose output includes extensive solver output.  Debugging Residual output will print the complete residual vector to the screen.";

  str_fileout="Select files you would like to generate in this run. \n Basic output gives density profiles, g(r) if requested, and principle output (adsorption and free energies). \nExtended output also gives external field files, interaction potential files, and polymer G functions. \nVerbose output gives additional mesh information, stencils, parallel computing informaiton, timings, and generates density profiles at each iteration. \nDebug Matrix prints the entire matrix to files (which are very large).  The code exits after the matrix is printed.";

  str_printtooutput="Select the bulk fluid output you would like to have printed in the file dft_output.dat.  Note that density and pressure output is not available when chemical potential continuation is being done.";
  str_adsout="Indicate how you would like the adsorption output to be printed in dft_output.dat";
  str_energyout="Indicate how you would like the free energy output to be printed in dft_output.dat";
  str_LperArea="Set to true to convert all output (energy, adsorption, force) to per unit area. \n If false, output will be energy (3D), energy/Length (2D), and energy/Area (1D).";
  str_reflect="Set to true if you want output (energy/adsorption/force) to use reflected images in calculations.";
  str_pmf="Set to true to produce output for surface-surface interactions.\n Needed for potential of mean force calculations.";

  /**************************************/
  /* Define validators for this section.*/
  /**************************************/
   RCP<StringValidator> ScreenOutValidator = rcp(
       new StringValidator(tuple<std::string>("No Screen Output","Errors Only to Screen", "Basic Screen Output",
                                              "Verbose Screen Output","Debugging Output: Residual")));

   RCP<StringValidator> FileOutValidator = rcp(
       new StringValidator(tuple<std::string>("Basic Files","Extended Files", "Debug Files","Debug Matrix Files")));

   RCP<StringValidator> StateOutputValidator = rcp( new StringValidator(tuple<std::string>(
     "Density[i]","Betamu[i]","kappa","Density[i], Betamu[i], kappa, and pressure","Density[all], Betamu[all], kappa, and pressure","No State Point Output")));

   RCP<StringValidator> StateOutputNoChargeValidator = rcp( new StringValidator(tuple<std::string>(
     "Density[i]","Betamu[i]","Density[i], Betamu[i], and pressure","Density[all], Betamu[all], and pressure","No State Point Output")));

   RCP<StringValidator> AdsEnergyOutputValidator = rcp( new StringValidator(tuple<std::string>(
       "total adsorption and energy","excess adsorption and energy", "excess and total adsorption and energy","adsorption/volume and energy/volume (bulk density and pressure)")));

   RCP<StringValidator> EnergyOutputValidator = rcp( new StringValidator(tuple<std::string>(
       "total free energy","excess surface free energy", "excess and total free energy","free energy/volume (pressure in bulk fluid)")));

   RCP<StringValidator> MeshOutputValidator = rcp( new StringValidator(tuple<std::string>("Position of all surfaces","Separations between surfaces")));


  /********************/
  /* Setup Parameters */
  /********************/

   Output_List->set("O1: Screen Output", "Basic Screen Output", str_screenout,ScreenOutValidator);
   Output_List->set("O2: Files Output", "Basic Files", str_fileout,FileOutValidator);

   if (Functional_List->get<string>("F3_CHARGE_Functional") !="No Charge or No Poisson" || 
       Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut/shift)" ||
       Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut only)" ){ 
         Output_List->set("O3.0: dft_output.dat: State Point Output", "Density[all], Betamu[all], kappa, and pressure", str_printtooutput,StateOutputValidator);
   }
   else{ Output_List->set("O3.0: dft_output.dat: State Point Output", "Density[all], Betamu[all], and pressure", str_printtooutput, StateOutputNoChargeValidator); }

   Output_List->set("O3.1: dft_output.dat: Adsorption and Energy Output", "excess and total adsorption and energy", str_adsout, AdsEnergyOutputValidator); 
   Output_List->set("O3.2: dft_output.dat: per unit area?", false, str_LperArea);
   Output_List->set("O3.3: dft_output.dat: Correct for reflections?", true, str_reflect);
   Output_List->set("O3.4: dft_output.dat: Mesh Output", "Separations between surfaces", "Select output type for mesh continuation.",MeshOutputValidator);

   Output_List->set("O4: Print g(r)?", false, "Set to true to produce a radial correlation function.");
   Output_List->set("O5: Print surface-surface interactions?", false, str_pmf);


   return;
}
/****************************************************************************************************************/
void dft_GUI_OutputParams_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Output_List)
{
   string str_screenout,str_fileout,str_printtooutput,str_adsout,str_energyout,str_LperArea,str_reflect,str_pmf;

  /*******************************/
  /* Define a few useful strings */
  /*******************************/
  str_screenout="Select how much output to the screen you would like to generate.\n Basic output includes thermodynamics and post processing results, but minimal output from solver.  \nVerbose output includes extensive solver output.  Debugging Residual output will print the complete residual vector to the screen.";

  str_fileout="Select files you would like to generate in this run. \n Basic output gives density profiles, g(r) if requested, and principle output (adsorption and free energies). \nExtended output also gives external field files, interaction potential files, and polymer G functions. \nVerbose output gives additional mesh information, stencils, parallel computing informaiton, timings, and generates density profiles at each iteration. \nDebug Matrix prints the entire matrix to files (which are very large).  The code exits after the matrix is printed.";

  str_printtooutput="Select the bulk fluid output you would like to have printed in the file dft_output.dat.  Note that density and pressure output is not available when chemical potential continuation is being done.";
  str_adsout="Indicate how you would like the adsorption output to be printed in dft_output.dat";
  str_energyout="Indicate how you would like the free energy output to be printed in dft_output.dat";
  str_LperArea="Set to true to convert all output (energy, adsorption, force) to per unit area. \n If false, output will be energy (3D), energy/Length (2D), and energy/Area (1D).";
  str_reflect="Set to true if you want output (energy/adsorption/force) to use reflected images in calculations.";
  str_pmf="Set to true to produce output for surface-surface interactions.\n Needed for potential of mean force calculations.";

  /**************************************/
  /* Define validators for this section.*/
  /**************************************/
   RCP<StringValidator> ScreenOutValidator = rcp(
       new StringValidator(tuple<std::string>("No Screen Output","Errors Only to Screen", "Basic Screen Output",
                                              "Verbose Screen Output","Debugging Output: Residual")));

   RCP<StringValidator> FileOutValidator = rcp(
       new StringValidator(tuple<std::string>("Basic Files","Extended Files", "Debug Files","Debug Matrix Files")));

   RCP<StringValidator> StateOutputValidator = rcp( new StringValidator(tuple<std::string>(
     "Density[i]","Betamu[i]","kappa","Density[i], Betamu[i], kappa, and pressure","Density[all], Betamu[all], kappa, and pressure","No State Point Output")));

   RCP<StringValidator> StateOutputNoChargeValidator = rcp( new StringValidator(tuple<std::string>(
     "Density[i]","Betamu[i]","Density[i], Betamu[i], and pressure","Density[all], Betamu[all], and pressure","No State Point Output")));

   RCP<StringValidator> AdsEnergyOutputValidator = rcp( new StringValidator(tuple<std::string>(
       "total adsorption and energy","excess adsorption and energy", "excess and total adsorption and energy","adsorption/volume and energy/volume (bulk density and pressure)")));

   RCP<StringValidator> EnergyOutputValidator = rcp( new StringValidator(tuple<std::string>(
       "total free energy","excess surface free energy", "excess and total free energy","free energy/volume (pressure in bulk fluid)")));

   RCP<StringValidator> MeshOutputValidator = rcp( new StringValidator(tuple<std::string>("Position of all surfaces","Separations between surfaces")));


   /****************************************/
   /* Set up parameters from Tramonto data */
   /****************************************/

   if (Iwrite==MINIMAL || Iwrite==DENSITIES){
      Output_List->set("O1: Screen Output", "Basic Screen Output", str_screenout,ScreenOutValidator);
      Output_List->set("O2: Files Output", "Basic Files", str_fileout,FileOutValidator);
   }
   else if (Iwrite==EXTENDED){
      Output_List->set("O1: Screen Output", "Basic Screen Output", str_screenout,ScreenOutValidator);
      Output_List->set("O2: Files Output", "Extended Files", str_fileout,FileOutValidator);
   }
   else if (Iwrite==VERBOSE){
      Output_List->set("O1: Screen Output", "Verbose Screen Output", str_screenout,ScreenOutValidator);
      Output_List->set("O2: Files Output", "Debug Files", str_fileout,FileOutValidator);
   }
   else if (Iwrite==NO_SCREEN){
     Output_List->set("O1: Screen Output", "No Screen Output", str_screenout,ScreenOutValidator);
     Output_List->set("O2: Files Output", "Basic Files", str_fileout,FileOutValidator);
   }
   else if (Iwrite==VERBOSE_MATRIX){
      Output_List->set("O1: Screen Output", "Debugging Output: Residual", str_screenout,ScreenOutValidator);
      Output_List->set("O2: Files Output", "Debug Matrix Files", str_fileout,FileOutValidator);
   }

   if (Print_rho_switch==SWITCH_NO_STATEOUT){
         if (Functional_List->get<string>("F3_CHARGE_Functional") !="No Charge or No Poisson" || 
            Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut/shift)" ||
            Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut only)" ){ 
              Output_List->set("O3.0: dft_output.dat: State Point Output", "No State Point Output", str_printtooutput,StateOutputValidator);
         }
         else Output_List->set("O3.0: dft_output.dat: State Point Output", "No State Point Output", str_printtooutput,StateOutputNoChargeValidator);
   }
   else{
      if (Print_rho_switch==SWITCH_ALLTYPES || Print_rho_switch==SWITCH_BULK_OUTPUT_ALL){
         if (Functional_List->get<string>("F3_CHARGE_Functional") !="No Charge or No Poisson" || 
            Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut/shift)" ||
            Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut only)" ){ 
              Output_List->set("O3.0: dft_output.dat: State Point Output", "Density[all], Betamu[all], kappa, and pressure", str_printtooutput,StateOutputValidator);
         }
         else Output_List->set("O3.0: dft_output.dat: State Point Output", "Density[all], Betamu[all], and pressure", str_printtooutput,StateOutputNoChargeValidator);

      }
      else if (Print_rho_switch==SWITCH_BULK_OUTPUT || Print_rho_switch==SWITCH_ALLTYPES_ICOMP){
         if (Functional_List->get<string>("F3_CHARGE_Functional") !="No Charge or No Poisson" || 
            Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut/shift)" ||
            Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut only)" ){ 
              Output_List->set("O3.0: dft_output.dat: State Point Output", "Density[i], Betamu[i], kappa, and pressure", str_printtooutput,StateOutputValidator);
         }
         else Output_List->set("O3.0: dft_output.dat: State Point Output", "Density[i], Betamu[i], and pressure", str_printtooutput,StateOutputNoChargeValidator);
      }
      else if (Print_rho_switch==SWITCH_RHO){
         if (Functional_List->get<string>("F3_CHARGE_Functional") !="No Charge or No Poisson" || 
            Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut/shift)" ||
            Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut only)" ){ 
              Output_List->set("O3.0: dft_output.dat: State Point Output", "Density[i]", str_printtooutput,StateOutputValidator);
         }
         else Output_List->set("O3.0: dft_output.dat: State Point Output", "Density[i]", str_printtooutput,StateOutputNoChargeValidator);
      }
      else if (Print_rho_switch==SWITCH_ION)
            Output_List->set("O3.0: dft_output.dat: State Point Output", "kappa", str_printtooutput,StateOutputValidator);
      else if (Print_rho_switch==SWITCH_MU){
         if (Functional_List->get<string>("F3_CHARGE_Functional") !="No Charge or No Poisson" || 
            Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut/shift)" ||
            Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut only)" ){ 
              Output_List->set("O3.0: dft_output.dat: State Point Output", "Betamu[i]", str_printtooutput,StateOutputValidator);
         }
         else Output_List->set("O3.0: dft_output.dat: State Point Output", "Betamu[i]", str_printtooutput,StateOutputNoChargeValidator);
      } 
   }

   if (Functional_List->get<string>("F3_CHARGE_Functional") !="No Charge or No Poisson" || 
       Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut/shift)" ||
       Fluid_List->get<string>("F4_PairPotType") == "Coulomb potential as mean field (cut only)" ){ 
         Output_List->set("O3.0: dft_output.dat: State Point Output", "Density[all], Betamu[all], kappa, and pressure", str_printtooutput,StateOutputValidator);
   }
   else{ Output_List->set("O3.0: dft_output.dat: State Point Output", "Density[all], Betamu[all], and pressure", str_printtooutput, StateOutputNoChargeValidator); }

   if (Print_rho_switch==SWITCH_BULK_OUTPUT_ALL || Print_rho_switch == SWITCH_BULK_OUTPUT)
      Output_List->set("O3.1: dft_output.dat: Adsorption and Energy Output", "adsorption/volume and energy/volume (bulk density and pressure)", str_adsout,AdsEnergyOutputValidator); 
   else Output_List->set("O3.1: dft_output.dat: Adsorption and Energy Output", "excess and total adsorption and energy", str_adsout,AdsEnergyOutputValidator); 
   

   if (Lper_area==TRUE) Output_List->set("O3.2: dft_output.dat: per unit area?", true, str_LperArea);
   else                 Output_List->set("O3.2: dft_output.dat: per unit area?", false, str_LperArea);

   if (Lcount_reflect==TRUE) Output_List->set("O3.3: dft_output.dat: Correct for reflections?", true, str_reflect);
   else Output_List->set("O3.3: dft_output.dat: Correct for reflections?", true, str_reflect);

   if (Print_mesh_switch==SWITCH_SURFACE_SEP) Output_List->set("O3.4: dft_output.dat: Mesh Output", "Separations between surfaces", "Select output type for mesh continuation.",MeshOutputValidator);
   else Output_List->set("O3.4: dft_output.dat: Mesh Output", "Position of all surfaces", "Select output type for mesh continuation.",MeshOutputValidator);

   if (Lprint_gofr==TRUE) Output_List->set("O4: Print g(r)?", true, "Set to true to produce a radial correlation function.");
   else Output_List->set("O4: Print g(r)?", false, "Set to true to produce a radial correlation function.");

   if (Lprint_pmf==TRUE) Output_List->set("O5: Print surface-surface interactions?", true, str_pmf);
   else Output_List->set("O5: Print surface-surface interactions?", false, str_pmf);

   return;
}
/****************************************************************************************************************/
void dft_GUI_OutputParams_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Mesh_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Continuation_List, 
                         Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Output_List)
{
    /***********************/
    /* set up dependencies */
    /***********************/
    RCP<StringVisualDependency> GofRprint_Dep = rcp(
         new StringVisualDependency( SurfaceInteraction_List->getEntryRCP("SI2.0: Vext selections"),Output_List->getEntryRCP("O4: Print g(r)?"),  
             tuple<std::string>("Vext=U(r) only (fixed atoms or colloids)")));

    RCP<StringVisualDependency> UWWprint_Dep = rcp(
         new StringVisualDependency( SurfaceInteraction_List->getEntryRCP("SI2.0: Vext selections"),Output_List->getEntryRCP("O5: Print surface-surface interactions?"),
             tuple<std::string>("Vext=U(r) only (fixed atoms or colloids)")));

/*      RCP<NumberVisualDependency<int> > LperAreaDep = rcp(
         new NumberVisualDependency<int>( Mesh_List->getEntryRCP("M1_Ndim"), 
                                          Output_List->getEntryRCP("O3.2: dft_output.dat: per unit area?"),true,
                                          Mesh_List->get<int>("M1_Ndim")-1));*/

    RCP<StringVisualDependency> MeshCont_Dep = rcp(
         new StringVisualDependency(Continuation_List->getEntryRCP("C1: Continuation Type"),Output_List->getEntryRCP("O3.4: dft_output.dat: Mesh Output"),
             tuple<std::string>("Mesh Continuation","Mesh Stepping with LOCA Binodal")));

    depSheet_Tramonto->addDependency(GofRprint_Dep);
    depSheet_Tramonto->addDependency(UWWprint_Dep);
/*      depSheet_Tramonto->addDependency(LperArea_Dep);*/
    depSheet_Tramonto->addDependency(MeshCont_Dep);

    return;
}
/****************************************************************************************************************/

