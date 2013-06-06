//using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_DensityStartupParams_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> DensProfile_List)
{
   string str_densityType, str_densityConstructor,str_fieldConstructor, str_extfieldType;
   string densFile1_default,densFile2_default;

   densFile1_default=(string)Runpath+"dft_dens.dat";
   densFile2_default=(string)Runpath+"dft_dens2.dat";

   /************************************/
   /* a few helpful string constants */
   /************************************/

   str_densityType="Select an option for the construction of a density profile.",
   str_fieldConstructor="Select an option for the construction of dependent fields (nonlocal densities, chain equations etc)";
   str_extfieldType="Select an option for setting up or reading in an external field.  May be summed from up to two files";

    /**************************************/
    /* Define validators for this section.*/
    /**************************************/
    RCP<StringValidator> GuessDensValidator = rcp(
        new StringValidator(tuple<std::string>("Constant Bulk Density","Rho_bulk*exp(-Vext/kT)", "Step function profile",
                          "Chopped profile to rho_bulk","Chopped profile to rho_step","Linear profile (for diffusion)")));

    RCP<StringValidator> GuessFieldsValidator = rcp(
        new StringValidator(tuple<std::string>("Bulk values for dependent fields","Compute fields based on density profile", 
              "Compute nonlocal densities only - other fields are bulk","Compute chain eq. and nonlocal densities - other fields are bulk")));

    RCP<StringValidator> RestartValidator = rcp( new StringValidator(tuple<std::string>(
      "Construct a new Density Profile","Restart from File","Restart with Step to constant",
      "Restart densities only (no other fields)", "Restart incomplete Ncomp","Restart with 1D profile (in 2D or 3D)")));

    RCP<StringValidator> RestartVextValidator = rcp( new StringValidator(tuple<std::string>(
      "Compute new Vext","Restart from File","Sum Two from Files", "Sum Two with constraints")));

   RCP<FileNameValidator> VextFileValidator = rcp(new FileNameValidator);
   RCP<FileNameValidator> VextFile2Validator = rcp(new FileNameValidator);

   RCP<FileNameValidator> DensityFileValidator = rcp(new FileNameValidator);
   RCP<FileNameValidator> DensityFile2Validator = rcp(new FileNameValidator);


   DensProfile_List->set("IG1: Initial Guess Type", "Construct a new Density Profile", str_densityType,RestartValidator);
   DensProfile_List->set("IG1.1: Density Restart File", densFile1_default, "Select the file that contains a density profile for restart", DensityFileValidator);
   DensProfile_List->set("IG1.2: 2nd Density RestartFile", densFile2_default, "Select a 2nd file that contains a density profile for binodal restart", DensityFile2Validator);
   DensProfile_List->set("IG1.3: Number of missing components", 0, "Set the number of components that are not included in the file.  They must have the highest density IDs");
   DensProfile_List->set("IG1.4: Rho max", 1000., "Set the maximum value allowed in densities that are read from a file." );

   DensProfile_List->set("IG2: Density Profile Constructor", "Constant Bulk Density", str_densityConstructor,GuessDensValidator);
   DensProfile_List->set("IG2.1: Nsteps", 2,"Set the number of steps in the constructed density profile");
   Array<int>Orient_step_Array(DensProfile_List->get<int>("IG2.1: Nsteps"),0);
   DensProfile_List->set("IG2.2: Orient_step[istep]", Orient_step_Array, "Set the orientation of each step in the profile.");

   Array<double>Xstart_step_Array(DensProfile_List->get<int>("IG2.1: Nsteps"),0.0);
   DensProfile_List->set("IG2.3: Xstart_step[istep]", Xstart_step_Array, "Set the starting position for each step in the profile. Note that the origin is in the center of the domain");

   Array<double>Xend_step_Array(DensProfile_List->get<int>("IG2.1: Nsteps"),0.0);
   DensProfile_List->set("IG2.4: Xend_step[istep]", Xend_step_Array, "Set the ending position for each step in the profile. Note that the origin is in the center of the domain");

   TwoDArray<double> RhoStep_Array(Fluid_List->get<int>("F1_Ncomp"),DensProfile_List->get<int>("IG2.1: Nsteps"),0.1);
   DensProfile_List->set("IG2.5 Rho_step[icomp][istep]", RhoStep_Array, "Defined densities for stepped profile.");

   DensProfile_List->set("IG3: Dependent Field Constructor", "Compute fields based on density profile", str_fieldConstructor,GuessFieldsValidator);
   DensProfile_List->set("IG4: External Field Constructor", "Compute new Vext", str_extfieldType,RestartVextValidator);

   DensProfile_List->set("IG4.1: External Field Filename", "dft_vext.dat", "Select the file that contains an external field for restart", VextFileValidator);
   DensProfile_List->set("IG4.2: 2nd External Field Filename", "", "Select a 2nd (static to continuation) file that contains an external field for restart", VextFile2Validator);

   return;
}
/*******************************************************************************************************************/
void dft_GUI_DensityStartupParams_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> DensProfile_List)
{
   int i,j;
   string str_densityType, str_densityConstructor,str_fieldConstructor, str_extfieldType;
   string densFile1_default,densFile2_default;

   densFile1_default=(string)Runpath+"dft_dens.dat";
   densFile2_default=(string)Runpath+"dft_dens2.dat";

   /************************************/
   /* a few helpful string constants */
   /************************************/

   str_densityType="Select an option for the construction of a density profile.",
   str_fieldConstructor="Select an option for the construction of dependent fields (nonlocal densities, chain equations etc)";
   str_extfieldType="Select an option for setting up or reading in an external field.  May be summed from up to two files";

    /**************************************/
    /* Define validators for this section.*/
    /**************************************/
    RCP<StringValidator> GuessDensValidator = rcp(
        new StringValidator(tuple<std::string>("Constant Bulk Density","Rho_bulk*exp(-Vext/kT)", "Step function profile",
                          "Chopped profile to rho_bulk","Chopped profile to rho_step","Linear profile (for diffusion)")));

    RCP<StringValidator> GuessFieldsValidator = rcp(
        new StringValidator(tuple<std::string>("Bulk values for dependent fields","Compute fields based on density profile", 
              "Compute nonlocal densities only - other fields are bulk","Compute chain eq. and nonlocal densities - other fields are bulk")));

    RCP<StringValidator> RestartValidator = rcp( new StringValidator(tuple<std::string>(
      "Construct a new Density Profile","Restart from File","Restart with Step to constant",
      "Restart densities only (no other fields)", "Restart incomplete Ncomp","Restart with 1D profile (in 2D or 3D)")));

    RCP<StringValidator> RestartVextValidator = rcp( new StringValidator(tuple<std::string>(
      "Compute new Vext","Restart from File","Sum Two from Files", "Sum Two with constraints")));

   RCP<FileNameValidator> VextFileValidator = rcp(new FileNameValidator);
   RCP<FileNameValidator> VextFile2Validator = rcp(new FileNameValidator);

   RCP<FileNameValidator> DensityFileValidator = rcp(new FileNameValidator);
   RCP<FileNameValidator> DensityFile2Validator = rcp(new FileNameValidator);

   if (Restart==NORESTART) DensProfile_List->set("IG1: Initial Guess Type", "Construct a new Density Profile", str_densityType,RestartValidator);
   else if (Restart==RESTART_BASIC) DensProfile_List->set("IG1: Initial Guess Type", "Restart from File", str_densityType,RestartValidator);
   else if (Restart==RESTART_STEP) DensProfile_List->set("IG1: Initial Guess Type", "Restart with Step to constant", str_densityType,RestartValidator);
   else if (Restart==RESTART_DENSONLY) DensProfile_List->set("IG1: Initial Guess Type", "Restart densities only (no other fields)", str_densityType,RestartValidator);
   else if (Restart==RESTART_FEWERCOMP) DensProfile_List->set("IG1: Initial Guess Type", "Restart incomplete Ncomp", str_densityType,RestartValidator);
   else if (Restart==RESTART_1DTOND) DensProfile_List->set("IG1: Initial Guess Type", "Restart with 1D profile (in 2D or 3D)", str_densityType,RestartValidator);

   DensProfile_List->set("IG1.1: Density Restart File",densFile1_default, "Select the file that contains an external field for restart", DensityFileValidator);
   DensProfile_List->set("IG1.2: 2nd Density RestartFile",densFile2_default, "Select a 2nd (static to continuation) file that contains an external field for restart", DensityFile2Validator);

   DensProfile_List->set("IG1.3: Number of missing components", Nmissing_densities, "Set the number of components that are not included in the file.  They must have the highest density IDs");
   DensProfile_List->set("IG1.4: Rho max", Rho_max, "Set the maximum value allowed in densities that are read from a file." );

   if (Iguess==CONST_RHO) DensProfile_List->set("IG2: Density Profile Constructor", "Constant Bulk Density", str_densityConstructor,GuessDensValidator);
   else if (Iguess==EXP_RHO) DensProfile_List->set("IG2: Density Profile Constructor", "Rho_bulk*exp(-Vext/kT)", str_densityConstructor,GuessDensValidator);
   else if (Iguess==STEP_PROFILE) DensProfile_List->set("IG2: Density Profile Constructor", "Step function profile", str_densityConstructor,GuessDensValidator);
   else if (Iguess==CHOP_RHO) DensProfile_List->set("IG2: Density Profile Constructor", "Chopped profile to rho_bulk", str_densityConstructor,GuessDensValidator);
   else if (Iguess==CHOP_RHO_STEP) DensProfile_List->set("IG2: Density Profile Constructor", "Chopped profile to rho_step", str_densityConstructor,GuessDensValidator);
   else if (Iguess==LINEAR) DensProfile_List->set("IG2: Density Profile Constructor", "Linear profile (for diffusion)", str_densityConstructor,GuessDensValidator);

   DensProfile_List->set("IG2.1: Nsteps", Nsteps,"Set the number of steps in the constructed density profile");

   Array<int>Orient_step_Array(Orientation_step,Orientation_step+Nsteps);
   DensProfile_List->set("IG2.2: Orient_step[istep]", Orient_step_Array, "Set the orientation of each step in the profile.");

   Array<double>Xstart_step_Array(Xstart_step,Xstart_step+Nsteps);
   DensProfile_List->set("IG2.3: Xstart_step[istep]", Xstart_step_Array, "Set the starting position for each step in the profile. Note that the origin is in the center of the domain");

   Array<double>Xend_step_Array(Xend_step,Xend_step+Nsteps);
   DensProfile_List->set("IG2.4: Xend_step[istep]", Xend_step_Array, "Set the ending position for each step in the profile. Note that the origin is in the center of the domain");

   TwoDArray<double> RhoStep_Array(Ncomp,Nsteps);
   for (i=0;i<Ncomp;i++) for(j=0;j<Nsteps;j++) RhoStep_Array(i,j)=Rho_step[i][j];
   DensProfile_List->set("IG2.5 Rho_step[icomp][istep]", RhoStep_Array, "Defined densities for stepped profile.");

   if (Iguess_fields==BULK) DensProfile_List->set("IG3: Dependent Field Constructor", "Bulk values for dependent fields", str_fieldConstructor,GuessFieldsValidator);
   else if (Iguess_fields==CALC_ALL_FIELDS) DensProfile_List->set("IG3: Dependent Field Constructor", "Compute fields based on density profile", str_fieldConstructor,GuessFieldsValidator);
   else if (Iguess_fields==CALC_RHOBAR_ONLY) DensProfile_List->set("IG3: Dependent Field Constructor", "Compute nonlocal densities only - other fields are bulk", str_fieldConstructor,GuessFieldsValidator);
   else if (Iguess_fields==CALC_RHOBAR_AND_G) DensProfile_List->set("IG3: Dependent Field Constructor", "Compute chain eq. and nonlocal densities - other fields are bulk", str_fieldConstructor,GuessFieldsValidator);

   if (Restart_Vext==READ_VEXT_FALSE) DensProfile_List->set("IG4: External Field Constructor", "Compute new Vext", str_extfieldType,RestartVextValidator);
   else if (Restart_Vext==READ_VEXT_TRUE) DensProfile_List->set("IG4: External Field Constructor", "Restart from File", str_extfieldType,RestartVextValidator);
   else if (Restart_Vext==READ_VEXT_SUMTWO) DensProfile_List->set("IG4: External Field Constructor", "Sum Two from Files", str_extfieldType,RestartVextValidator);
   else if (Restart_Vext==READ_VEXT_STATIC) DensProfile_List->set("IG4: External Field Constructor", "Sum Two with constraints", str_extfieldType,RestartVextValidator);

   if (Vext_filename!=NULL) 
        DensProfile_List->set("IG4.1: External Field Filename",(string)Vext_filename, "Select the file that contains an external field for restart", VextFileValidator);
   else DensProfile_List->set("IG4.1: External Field Filename","", "Select the file that contains an external field for restart", VextFileValidator);
   if (Vext_filename2!=NULL) 
        DensProfile_List->set("IG4.2: 2nd External Field Filename",(string)Vext_filename2, "Select a 2nd (static to continuation) file that contains an external field for restart", VextFile2Validator);
    else DensProfile_List->set("IG4.2: 2nd External Field Filename","", "Select a 2nd (static to continuation) file that contains an external field for restart", VextFile2Validator);

    return;
}
/*******************************************************************************************************************/
void dft_GUI_DensityStartupParams_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> Continuation_List,
                         Teuchos::RCP<Teuchos::ParameterList> DensProfile_List)
{


    RCP<StringVisualDependency> Nmissing_Dep = rcp(
         new StringVisualDependency( DensProfile_List->getEntryRCP("IG1: Initial Guess Type"),DensProfile_List->getEntryRCP("IG1.3: Number of missing components"),  
             tuple<std::string>("Restart incomplete Ncomp")));

    Dependency::ParameterEntryList Restart_Deps;
    Restart_Deps.insert(DensProfile_List->getEntryRCP("IG1.4: Rho max"));
    Restart_Deps.insert(DensProfile_List->getEntryRCP("IG1.1: Density Restart File"));
    RCP<StringVisualDependency> Restart_Dep = rcp(
         new StringVisualDependency( DensProfile_List->getEntryRCP("IG1: Initial Guess Type"),Restart_Deps,"Construct a new Density Profile",false));

    RCP<StringCondition> Cont_Cond=rcp(new StringCondition(Continuation_List->getEntryRCP("C1: Continuation Type"),
                          tuple<std::string>("LOCA: Binodal Continuation","Mesh Stepping with LOCA Binodal","LOCA: Spinodal Continuation")));
    RCP<StringCondition> Start_Cond=rcp(new StringCondition(DensProfile_List->getEntryRCP("IG1: Initial Guess Type"),
                           tuple<std::string>("Restart from File","Restart with Step to constant",
                                              "Restart densities only (no other fields)","Restart incomplete Ncomp",
                                              "Restart with 1D profile (in 2D or 3D)")));
    Condition::ConstConditionList DensProf2_ConList=tuple<RCP<const Condition> >(Cont_Cond,Start_Cond);
    RCP<AndCondition>BinodalRead_Con=rcp(new AndCondition(DensProf2_ConList));
    RCP<ConditionVisualDependency> BinodalRead_Dep = rcp(
         new ConditionVisualDependency(BinodalRead_Con,DensProfile_List->getEntryRCP("IG1.2: 2nd Density RestartFile"),true));

    RCP<StringCondition> ProfStCond1=rcp(new StringCondition(DensProfile_List->getEntryRCP("IG1: Initial Guess Type"),"Construct a new Density Profile"));
    RCP<StringCondition> ProfStCond2=rcp(new StringCondition(DensProfile_List->getEntryRCP("IG1: Initial Guess Type"),"Restart incomplete Ncomp"));
    Condition::ConstConditionList ProfStep_ConList=tuple<RCP<const Condition> >(ProfStCond1,ProfStCond2);
    RCP<OrCondition>ProfConstruct_orCon=rcp(new OrCondition(ProfStep_ConList));

    RCP<ConditionVisualDependency> DensCreate_Dep = rcp(
         new ConditionVisualDependency(ProfConstruct_orCon,DensProfile_List->getEntryRCP("IG2: Density Profile Constructor"),true));

    RCP<StringCondition> StepOptions_Cond1=rcp(new StringCondition(DensProfile_List->getEntryRCP("IG1: Initial Guess Type"),"Restart with Step to constant"));

    RCP<StringCondition> StepOptions_Cond2=rcp(new StringCondition(DensProfile_List->getEntryRCP("IG2: Density Profile Constructor"),
                               tuple<std::string>("Step function profile","Chopped profile to rho_bulk","Chopped profile to rho_step")));

    Condition::ConstConditionList StepOptions_ConList=tuple<RCP<const Condition> >(StepOptions_Cond1,StepOptions_Cond2);
    RCP<OrCondition>StepOptions_orCon=rcp(new OrCondition(StepOptions_ConList));

    Dependency::ParameterEntryList StepProfiles_Deps;
    StepProfiles_Deps.insert(DensProfile_List->getEntryRCP("IG2.1: Nsteps"));
    StepProfiles_Deps.insert(DensProfile_List->getEntryRCP("IG2.2: Orient_step[istep]"));
    StepProfiles_Deps.insert(DensProfile_List->getEntryRCP("IG2.3: Xstart_step[istep]"));
    StepProfiles_Deps.insert(DensProfile_List->getEntryRCP("IG2.4: Xend_step[istep]"));
    StepProfiles_Deps.insert(DensProfile_List->getEntryRCP("IG2.5 Rho_step[icomp][istep]"));

    RCP<ConditionVisualDependency> StepOptions_Dep = rcp( new ConditionVisualDependency(StepOptions_orCon,StepProfiles_Deps,true));

    RCP<StringVisualDependency> VextFile1_Dep = rcp(
         new StringVisualDependency( DensProfile_List->getEntryRCP("IG4: External Field Constructor"),DensProfile_List->getEntryRCP("IG4.1: External Field Filename"),  
             "Compute new Vext",false));

    RCP<StringVisualDependency> VextFile2_Dep = rcp(
         new StringVisualDependency( DensProfile_List->getEntryRCP("IG4: External Field Constructor"),DensProfile_List->getEntryRCP("IG4.2: 2nd External Field Filename"),  
             tuple<std::string>("Sum Two from Files","Sum Two with constraints")));


      /* put dependencies into dependency sheet */
    depSheet_Tramonto->addDependency(Nmissing_Dep);
    depSheet_Tramonto->addDependency(Restart_Dep);
    depSheet_Tramonto->addDependency(BinodalRead_Dep);
    depSheet_Tramonto->addDependency(DensCreate_Dep);
    depSheet_Tramonto->addDependency(StepOptions_Dep);
    depSheet_Tramonto->addDependency(VextFile1_Dep);
    depSheet_Tramonto->addDependency(VextFile2_Dep);

    return;
}
/*******************************************************************************************************************/
