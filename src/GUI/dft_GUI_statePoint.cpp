//using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_StatePoint_set_defaults( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                  Teuchos::RCP<Teuchos::ParameterList> Polymer_List, 
                  Teuchos::RCP<Teuchos::ParameterList> StatePoint_List,
                  Teuchos::RCP<Teuchos::ParameterList> Diffusion_List,
                  Teuchos::RCP<Teuchos::ParameterList> ChargedFluid_List)
{

 /**************************************/
 /* Define validators for this section.*/
 /**************************************/

  RCP<EnhancedNumberValidator<int> > DimValidator = rcp(new EnhancedNumberValidator<int>(0,2,1));
  RCP<StringValidator> DielecType_Validator = rcp(
         new StringValidator(tuple<std::string>("Uniform Dielectric Constant",
                   "2-Dielec Const: fluid and wall regions",
                   "3-Dielec Const: bulk fluid, fluid near wall, wall regions")));

  RCP<StringValidator> ElecPotType_Validator = rcp(
         new StringValidator(tuple<std::string>("Enter Dimensionless Potentials (phi e/kT)",
                    "Enter electrostatic potentials in mV units")));
 
  /*********************/
  /* set up parameters */
  /*********************/

  StatePoint_List->set("BF1: Dimensionless Density Entry?", true, "Dimensionless densities are given in units of rho*sigma_ref^3.\n Set to false to enter in any other units (molar,g/cc etc).\n For dimensioned entry you must provide a conversion factor, Fac, where rho*sigma*3=(your densities)/Fac.");
  StatePoint_List->set("BF2: Density Conversion Factor",1.0,"Provide a conversion factor, Fac, where rho*sigma*3=(your densities)/Fac.");

  Array<double> RhoBulk0_Array( (Fluid_List->get<int>("F1_Ncomp")),0.0);
  Array<double> RhoBulk0POL_Array( (Polymer_List->get<int>("P1: Npoly_comp")),0.0);
  StatePoint_List->set("BF3: Rho_b_0[icomp]", RhoBulk0_Array, "Set the bulk densities (homogeneous system), or the densities on the negative boundary (inhomogenous system)");
  StatePoint_List->set("BF3: Rho_b_0[ipol_comp]", RhoBulk0POL_Array, "Set the bulk densities (homogeneous system), or the densities on the negative boundary (inhomogenous system)");

  Array<double> RhoBulk1_Array( (Fluid_List->get<int>("F1_Ncomp")),0.0);
  Array<double> RhoBulk1POL_Array( (Polymer_List->get<int>("P1: Npoly_comp")),0.0);
  StatePoint_List->set("BF4: Rho_b_1[icomp]", RhoBulk1_Array, "Set the densities on the positive boundary  for an inhomogenous system");
  StatePoint_List->set("BF4: Rho_b_1[ipol_comp]", RhoBulk1POL_Array, "Set the densities on the positive boundary  for an inhomogenous system");

  StatePoint_List->set("BF5: Direction of Gradient", 0, "Indicate direction where inhomogeneous boundaries should be applied (0=x,1=y,2=z)",DimValidator);
  StatePoint_List->set("BF6: Constrained Interface?", false, "Set to true to set rho[0]=0.5(rho[0]_left+rho[0]_right)\n at the midpoint of the domain. This approach can help to keep the interface from moving in a free interface calculation.");

  StatePoint_List->set("BF7: X_const", 0.0, "Set the distance from the computational boundaries\n in the direction of Grad_dim\n where the chemical potential is to be held constant.\n Facilitates comparison with GCMD calculations.");


       /* Optional Diffusion parameters */
  Array<double> DiffCoeff_Array( (Fluid_List->get<int>("F1_Ncomp")),0.0);
  Array<double> DiffCoeffPOL_Array( (Polymer_List->get<int>("P1: Npoly_comp")),0.0);
  Diffusion_List->set("D1: Diff_Coeff[icomp]", DiffCoeff_Array, "Set the diffusion coefficientf sof each component.");
  Diffusion_List->set("D1: Diff_Coeff[ipol_comp]", DiffCoeffPOL_Array, "Set the diffusion coefficientf sof each component.");
  Diffusion_List->set("D3: Velocity", 0.0, "Set a bulk velocity term in diffusion calculation.");


       /* Bulk Charged Fluid Parameters */
  ChargedFluid_List->set("CF1: Type Dielectric Constant(s)","Uniform Dielectric Constant","Indicate type of dielectric constants to be entered for the problem.",DielecType_Validator);
  ChargedFluid_List->set("CF2: Entry of Relative Dielectric Constant(s)?",true,"true indicates that dielectric constants will be entered as D/D_reference\n (for example, D_ref may be the dielectric constant of the bulk solvent).\n If set to false, you will be asked for provide a reference value, D_reference for computation of relative dielectric constants.");
  ChargedFluid_List->set("CF3: Reference Dielectric Constant",78.5,"Enter a reference dielectric constant for calculation of relative values D/D_reference.\n Note that 78.5 is a typical value used for water");

  ChargedFluid_List->set("CF4.0: Dielec Const Bulk Fluid",1.0,"Enter the dielectric constant of the bulk fluid.");
  ChargedFluid_List->set("CF4.1: Dielec Const Near Wall Fluid",1.0,"Enter the dielectric constant of the fluid in the pore or near the wall.");
  ChargedFluid_List->set("CF4.2: Size of Near Wall region",1.0,"Enter the distance from the wall to be considered in the Near Wall region.");
  ChargedFluid_List->set("CF5.0: Sigma for Plasma Parameter (Angstroms)",4.25,"Enter the length scale you would like to use for the plasma parameter");
  ChargedFluid_List->set("CF5.1: Temperature for Plasma Parameter (Kelvin)",298.,"Enter the temperature you would like to use for calculation of the plasma parameter");
  ChargedFluid_List->set("CF5.2: Dielectric const for Plasma Parameter",KAPPA_H2O,"Enter the dielectric constant you would like to use for calculation of the plasma parameter");
  ChargedFluid_List->set("CF6.0: Electrostatic Potential(s) entry type","Enter Dimensionless Potentials (phi e/kT)","select how bulk fluid electrostatic potentials will be entered.\n  Note that setting phi=0 in a uniform bulk fluid is customary in electrostatics, but not for diffusive elctrochemical systems.",ElecPotType_Validator);
  ChargedFluid_List->set("CF6.1: Elec_pot_0", 0.0, "Set a bulk electrostatic potential (homogeneous boundary case), or \n set the bulk electrostatic potential on the negative boundary (inhomogeneous boundary or diffusion).");
  ChargedFluid_List->set("CF6.2: Elec_pot_1", 0.0, "Set the bulk electrostatic potential on the positive boundary (inhomogeneous boundary or diffusion).");

  return;
}
/****************************************************************************************************************/
void dft_GUI_StatePoint_set_OldFormat( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                  Teuchos::RCP<Teuchos::ParameterList> Polymer_List, 
                  Teuchos::RCP<Teuchos::ParameterList> StatePoint_List,
                  Teuchos::RCP<Teuchos::ParameterList> Diffusion_List,
                  Teuchos::RCP<Teuchos::ParameterList> ChargedFluid_List)
{
  string tmp_string;

  /**************************************/
  /* Define validators for this section.*/
  /**************************************/

   RCP<EnhancedNumberValidator<int> > DimValidator = rcp(new EnhancedNumberValidator<int>(0,2,1));
   RCP<StringValidator> DielecType_Validator = rcp(
         new StringValidator(tuple<std::string>("Uniform Dielectric Constant",
                   "2-Dielec Const: fluid and wall regions",
                   "3-Dielec Const: bulk fluid, fluid near wall, wall regions")));

   RCP<StringValidator> ElecPotType_Validator = rcp(
         new StringValidator(tuple<std::string>("Enter Dimensionless Potentials (phi e/kT)",
                    "Enter electrostatic potentials in mV units")));
 
   if (Density_ref<0.0){
        StatePoint_List->set("BF1: Dimensionless Density Entry?", true, "Dimensionless densities are given in units of rho*sigma_ref^3.\n Set to false to enter in any other units (molar,g/cc etc).\n For dimensioned entry you must provide a conversion factor, Fac, where rho*sigma*3=(your densities)/Fac.");
        StatePoint_List->set("BF2: Density Conversion Factor",1.0,"Provide a conversion factor, Fac, where rho*sigma*3=(your densities)/Fac.");
   }
   else{
        StatePoint_List->set("BF1: Dimensionless Density Entry?", false, "Dimensionless densities are given in units of rho*sigma_ref^3.\n Set to false to enter in any other units (molar,g/cc etc).\n For dimensioned entry you must provide a conversion factor, Fac, where rho*sigma*3=(your densities)/Fac.");
        StatePoint_List->set("BF2: Density Conversion Factor",Density_ref,"Provide a conversion factor, Fac, where rho*sigma*3=(your densities)/Fac.");
        
   }
   if (Type_interface != UNIFORM_INTERFACE){
      Array<double> RhoBulk0_Array(Rho_b_LBB,Rho_b_LBB+Ncomp);
      Array<double> RhoBulk0POL_Array(Rho_b_LBB,Rho_b_LBB+Npol_comp);
      StatePoint_List->set("BF3: Rho_b_0[icomp]", RhoBulk0_Array, "Set the bulk densities (homogeneous system), or the densities on the negative boundary (inhomogenous system)");
      StatePoint_List->set("BF3: Rho_b_0[ipol_comp]", RhoBulk0POL_Array, "Set the bulk densities (homogeneous system), or the densities on the negative boundary (inhomogenous system)");

      Array<double> RhoBulk1_Array(Rho_b_RTF,Rho_b_RTF+Ncomp);
      Array<double> RhoBulk1POL_Array(Rho_b_RTF,Rho_b_RTF+Npol_comp);
      StatePoint_List->set("BF4: Rho_b_1[icomp]", RhoBulk1_Array, "Set the densities on the positive boundary  for an inhomogenous system");
      StatePoint_List->set("BF4: Rho_b_1[ipol_comp]", RhoBulk1POL_Array, "Set the densities on the positive boundary  for an inhomogenous system");

      StatePoint_List->set("BF5: Direction of Gradient", Grad_dim, "Indicate direction where inhomogeneous boundaries should be applied (0=x,1=y,2=z)",DimValidator);

      if (Lconstrain_interface==TRUE)
            StatePoint_List->set("BF6: Constrained Interface?", true, "Set to true to set rho[0]=0.5(rho[0]_left+rho[0]_right)\n at the midpoint of the domain. This approach can help to keep the interface from moving in a free interface calculation.");
      else  StatePoint_List->set("BF6: Constrained Interface?", false, "Set to true to set rho[0]=0.5(rho[0]_left+rho[0]_right)\n at the midpoint of the domain. This approach can help to keep the interface from moving in a free interface calculation.");

   }
   else{
      Array<double> RhoBulk0_Array(Rho_b,Rho_b+Ncomp);
      Array<double> RhoBulk0POL_Array(Rho_b,Rho_b+Npol_comp);
      StatePoint_List->set("BF3: Rho_b_0[icomp]", RhoBulk0_Array, "Set the bulk densities (homogeneous system), or the densities on the negative boundary (inhomogenous system)");
      StatePoint_List->set("BF3: Rho_b_0[ipol_comp]", RhoBulk0POL_Array, "Set the bulk densities (homogeneous system), or the densities on the negative boundary (inhomogenous system)");

      Array<double> RhoBulk1_Array(Rho_b,Rho_b+Ncomp);
      Array<double> RhoBulk1POL_Array(Rho_b,Rho_b+Npol_comp);
      StatePoint_List->set("BF4: Rho_b_1[icomp]", RhoBulk1_Array, "Set the densities on the positive boundary  for an inhomogenous system");
      StatePoint_List->set("BF4: Rho_b_1[ipol_comp]", RhoBulk1POL_Array, "Set the densities on the positive boundary  for an inhomogenous system");

      StatePoint_List->set("BF5: Direction of Gradient", 0, "Indicate direction where inhomogeneous boundaries should be applied (0=x,1=y,2=z)",DimValidator);

      StatePoint_List->set("BF6: Constrained Interface?", false, "Set to true to set rho[0]=0.5(rho[0]_left+rho[0]_right)\n at the midpoint of the domain. This approach can help to keep the interface from moving in a free interface calculation.");

   }
   StatePoint_List->set("BF7: X_const", X_const_mu, "Set the distance from the computational boundaries\n in the direction of Grad_dim\n where the chemical potential is to be held constant.\n Facilitates comparison with GCMD calculations.");

       /* Optional Diffusion parameters */
   if (Type_interface != UNIFORM_INTERFACE){
        if (Type_interface != UNIFORM_INTERFACE){
           Array<double> DiffCoeff_Array(D_coef,D_coef+Ncomp);
           Array<double> DiffCoeffPOL_Array(D_coef,D_coef+Npol_comp);
           Diffusion_List->set("D1: Diff_Coeff[icomp]", DiffCoeff_Array, "Set the diffusion coefficientf sof each component. Note this is only used to post-process flux.  It is not neede for calculation of profiles.");
           Diffusion_List->set("D1: Diff_Coeff[ipol_comp]", DiffCoeffPOL_Array, "Set the diffusion coefficientf sof each component. Note this is only used to post-process flux.  It is not neede for calculation of profiles.");
           
           Diffusion_List->set("D3: Velocity", Velocity, "Optional diffusion parameter: Set a bulk velocity term in diffusion calculation.");
       }
   }

       /* Bulk Charged Fluid Parameters */
   if (Type_dielec==DIELEC_CONST) tmp_string="Uniform Dielectric Constant";
   else if (Type_dielec==DIELEC_WF) tmp_string="2-Dielec Const: fluid and wall regions";
   else if (Type_dielec==DIELEC_WF_PORE) tmp_string="3-Dielec Const: bulk fluid, fluid near wall, wall regions";
   ChargedFluid_List->set("CF1: Type Dielectric Constant(s)",tmp_string,"Indicate type of dielectric constants to be entered for the problem.",DielecType_Validator);

   if (Dielec_ref<=0.0) ChargedFluid_List->set("CF2: Entry of Relative Dielectric Constant(s)?",true,"true indicates that dielectric constants will be entered as D/D_reference\n (for example, D_ref may be the dielectric constant of the bulk solvent).\n If set to false, you will be asked for provide a reference value, D_reference for computation of relative dielectric constants.");
   else ChargedFluid_List->set("CF2: Entry of Relative Dielectric Constant(s)?",false,"true indicates that dielectric constants will be entered as D/D_reference\n (for example, D_ref may be the dielectric constant of the bulk solvent).\n If set to false, you will be asked for provide a reference value, D_reference for computation of relative dielectric constants.");

   ChargedFluid_List->set("CF3: Reference Dielectric Constant",Dielec_ref,"Enter a reference dielectric constant for calculation of relative values D/D_reference.\n Note that 78.5 is a typical value used for water");
   ChargedFluid_List->set("CF4.0: Dielec Const Bulk Fluid",Dielec_bulk,"Enter the dielectric constant of the bulk fluid.");
   ChargedFluid_List->set("CF4.1: Dielec Const Near Wall Fluid",Dielec_pore,"Enter the dielectric constant of the fluid in the pore or near the wall.");
   ChargedFluid_List->set("CF4.2: Size of Near Wall region",Dielec_X,"Enter the distance from the wall to be considered in the Near Wall region.");
   ChargedFluid_List->set("CF5.0: Sigma for Plasma Parameter (Angstroms)",Sigma_Angstroms_plasma,"Enter the length scale you would like to use for calculation of the plasma parameter");
   ChargedFluid_List->set("CF5.1: Temperature for Plasma Parameter (Kelvin)",Temp_K_plasma,"Enter the temperature you would like to use for calculation of the plasma parameter");
   ChargedFluid_List->set("CF5.2: Dielectric const for Plasma Parameter",DielecConst_plasma,"Enter the dielectric constant you would like to use for calculation of the plasma parameter");

   ChargedFluid_List->set("CF6.0: Electrostatic Potential(s) entry type","Enter Dimensionless Potentials (phi e/kT)","select how bulk fluid electrostatic potentials will be entered.\n  Note that setting phi=0 in a uniform bulk fluid is customary in electrostatics, but not for diffusive elctrochemical systems.",ElecPotType_Validator);

   if (Type_interface != UNIFORM_INTERFACE){ 
      ChargedFluid_List->set("CF6.1: Elec_pot_0", Elec_pot_LBB, "Set a bulk electrostatic potential (homogeneous boundary case), or \n set the bulk electrostatic potential on the negative boundary (inhomogeneous boundary or diffusion).");
      ChargedFluid_List->set("CF6.2: Elec_pot_1", Elec_pot_RTF, "Set the bulk electrostatic potential on the positive boundary (inhomogeneous boundary or diffusion).");
   }  
   else{
      ChargedFluid_List->set("CF6.1: Elec_pot_0", 0.0, "Set a bulk electrostatic potential (homogeneous boundary case), or \n set the bulk electrostatic potential on the negative boundary (inhomogeneous boundary or diffusion).");
      ChargedFluid_List->set("CF6.2: Elec_pot_1", 0.0, "Set the bulk electrostatic potential on the positive boundary (inhomogeneous boundary or diffusion).");
   }

   return;
}
/****************************************************************************************************************/
void dft_GUI_StatePoint_dependencies( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                  Teuchos::RCP<Teuchos::ParameterList> Polymer_List, 
                  Teuchos::RCP<Teuchos::ParameterList> StatePoint_List,
                  Teuchos::RCP<Teuchos::ParameterList> Diffusion_List,
                  Teuchos::RCP<Teuchos::ParameterList> ChargedFluid_List)
{
   /************************/
   /* set up dependencies */
   /************************/

    RCP<StringCondition> IpolcompCon = rcp(
           new StringCondition(Functional_List->getEntryRCP("F4_POLYMER_Functional"), 
                              tuple<std::string>("Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT",
                              "Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)",
                              "Polymer_JDC_iSAFT(comp)")));

    RCP<StringCondition> IcompCon = rcp(
           new StringCondition(Functional_List->getEntryRCP("F4_POLYMER_Functional"), "No Bonded Fluids"));

    RCP<StringCondition> InterfaceCon = rcp(
           new StringCondition(Functional_List->getEntryRCP("F0_Type_of_Calculation"), 
                              tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)",
                                   "Steady State Diffusion (inhomogeneous boundaries)")));

   RCP<StringCondition> DiffusionCon = rcp(
           new StringCondition(Functional_List->getEntryRCP("F0_Type_of_Calculation"), 
                                   "Steady State Diffusion (inhomogeneous boundaries)"));

    Condition::ConstConditionList Icomp_conList=tuple<RCP<const Condition> >(InterfaceCon,IcompCon);
    RCP<AndCondition> Rhoicomp_Con = rcp(new AndCondition(Icomp_conList));
      RCP<ConditionVisualDependency> IcompVis_Dep = rcp(
          new ConditionVisualDependency(Rhoicomp_Con, StatePoint_List->getEntryRCP("BF4: Rho_b_1[icomp]")));

    Condition::ConstConditionList Ipolcomp_conList=tuple<RCP<const Condition> >(InterfaceCon,IpolcompCon);
    RCP<AndCondition> Rhoipolcomp_Con = rcp(new AndCondition(Ipolcomp_conList));
      RCP<ConditionVisualDependency> IpolcompVis_Dep = rcp(
          new ConditionVisualDependency(Rhoipolcomp_Con, StatePoint_List->getEntryRCP("BF4: Rho_b_1[ipol_comp]")));


    Condition::ConstConditionList IpolcompDiff_conList=tuple<RCP<const Condition> >(DiffusionCon,IpolcompCon);
    RCP<AndCondition> Dcoefipolcomp_Con = rcp(new AndCondition(IpolcompDiff_conList));
      RCP<ConditionVisualDependency> DcoefIpolcompVis_Dep = rcp(
          new ConditionVisualDependency(Dcoefipolcomp_Con, Diffusion_List->getEntryRCP("D1: Diff_Coeff[ipol_comp]"))); 


    Condition::ConstConditionList IcompDiff_conList=tuple<RCP<const Condition> >(DiffusionCon,IcompCon);
    RCP<AndCondition> Dcoeficomp_Con = rcp(new AndCondition(IcompDiff_conList));
      RCP<ConditionVisualDependency> DcoefIcompVis_Dep = rcp(
          new ConditionVisualDependency(Dcoeficomp_Con, Diffusion_List->getEntryRCP("D1: Diff_Coeff[icomp]")));

   
   RCP<StringVisualDependency> GradDim_Dep = rcp(
       new StringVisualDependency( Functional_List->getEntryRCP("F0_Type_of_Calculation"), 
                                   StatePoint_List->getEntryRCP("BF5: Direction of Gradient"), 
                                   tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)",
                                   "Steady State Diffusion (inhomogeneous boundaries)")));

   RCP<StringVisualDependency> Diffusion_Dep = rcp(
       new StringVisualDependency( Functional_List->getEntryRCP("F0_Type_of_Calculation"),Diffusion_List->getEntryRCP("D3: Velocity"),
           tuple<std::string>("Steady State Diffusion (inhomogeneous boundaries)")));

   RCP<BoolVisualDependency> DensityUnit_Dep = rcp( 
        new BoolVisualDependency( StatePoint_List->getEntryRCP("BF1: Dimensionless Density Entry?"),
                                  StatePoint_List->getEntryRCP("BF2: Density Conversion Factor"),false));

   RCP<StringVisualDependency> Lconstrain_Dep = rcp(
       new StringVisualDependency( Functional_List->getEntryRCP("F0_Type_of_Calculation"),
                                   StatePoint_List->getEntryRCP("BF6: Constrained Interface?"),
                                   tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)"))); 

   RCP<StringVisualDependency> Rhoicomp_Dep = rcp(
       new StringVisualDependency( Functional_List->getEntryRCP("F4_POLYMER_Functional"),
                                   StatePoint_List->getEntryRCP("BF3: Rho_b_0[icomp]"),
                                   tuple<std::string>("No Bonded Fluids")),true); 

   RCP<StringVisualDependency> Rhoipolcomp_Dep = rcp(
       new StringVisualDependency( Functional_List->getEntryRCP("F4_POLYMER_Functional"),
                                   StatePoint_List->getEntryRCP("BF3: Rho_b_0[ipol_comp]"),
                                   tuple<std::string>("Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT",
                                   "Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)",
                                   "Polymer_JDC_iSAFT(comp)")),true); 


   Dependency::ParameterEntryList DensityLengthDependents;
   DensityLengthDependents.insert(StatePoint_List->getEntryRCP("BF3: Rho_b_0[icomp]"));
   DensityLengthDependents.insert(StatePoint_List->getEntryRCP("BF4: Rho_b_1[icomp]"));
   DensityLengthDependents.insert(Diffusion_List->getEntryRCP("D1: Diff_Coeff[icomp]"));
   RCP<NumberArrayLengthDependency<int,double> > DensityLength_Dep = rcp(
           new NumberArrayLengthDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"), DensityLengthDependents));

   Dependency::ParameterEntryList DensityLengthDependents_POL;
   DensityLengthDependents_POL.insert(StatePoint_List->getEntryRCP("BF3: Rho_b_0[ipol_comp]"));
   DensityLengthDependents_POL.insert(StatePoint_List->getEntryRCP("BF4: Rho_b_1[ipol_comp]"));
   DensityLengthDependents_POL.insert(Diffusion_List->getEntryRCP("D1: Diff_Coeff[ipol_comp]"));
   RCP<NumberArrayLengthDependency<int,double> > DensityLength_DepPOL = rcp(
           new NumberArrayLengthDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"), DensityLengthDependents_POL));


   Dependency::ParameterEntryList ChargeParamsDependents;
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF1: Type Dielectric Constant(s)"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF2: Entry of Relative Dielectric Constant(s)?"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF3: Reference Dielectric Constant"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF4.0: Dielec Const Bulk Fluid"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF5.0: Sigma for Plasma Parameter (Angstroms)"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF5.1: Temperature for Plasma Parameter (Kelvin)"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF5.2: Dielectric const for Plasma Parameter"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF6.0: Electrostatic Potential(s) entry type"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF6.1: Elec_pot_0"));

   RCP<StringVisualDependency> ChargeParams_Dep = rcp(
       new StringVisualDependency( Functional_List->getEntryRCP("F3_CHARGE_Functional"),ChargeParamsDependents,
           tuple<std::string>("Charge_Mean_Field","Charge_with_DeltaC_RPM","Charge_with_DeltaC_General","Charge(MF)_with_Polarization")));

   RCP<BoolVisualDependency> DielecConstUnit_Dep = rcp( 
        new BoolVisualDependency( ChargedFluid_List->getEntryRCP("CF2: Entry of Relative Dielectric Constant(s)?"),
                                  ChargedFluid_List->getEntryRCP("CF3: Reference Dielectric Constant"),false)); 

   Dependency::ParameterEntryList DielecPoreDependents;
   DielecPoreDependents.insert(ChargedFluid_List->getEntryRCP("CF4.1: Dielec Const Near Wall Fluid"));
   DielecPoreDependents.insert(ChargedFluid_List->getEntryRCP("CF4.2: Size of Near Wall region"));
   RCP<StringVisualDependency> DielecPore_Dep = rcp(
       new StringVisualDependency( ChargedFluid_List->getEntryRCP("CF1: Type Dielectric Constant(s)"),DielecPoreDependents,
           tuple<std::string>("3-Dielec Const: bulk fluid, fluid near wall, wall regions"))); 

       /* set up two dependees for the electrostatic condition */
   RCP<StringCondition> ChargeCon1 = rcp(
        new StringCondition(Functional_List->getEntryRCP("F3_CHARGE_Functional"),  
                            tuple<std::string>("Charge_Mean_Field","Charge_with_DeltaC_RPM",
                            "Charge_with_DeltaC_General","Charge(MF)_with_Polarization")));

   RCP<StringCondition> ChargeCon2 = rcp(
        new StringCondition(Functional_List->getEntryRCP("F0_Type_of_Calculation"),  
           tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)",
                              "Steady State Diffusion (inhomogeneous boundaries)")));

   Condition::ConstConditionList Charge_conList = tuple<RCP<const Condition> >(ChargeCon1, ChargeCon2);
   RCP<AndCondition> Charge_andCon = rcp(new AndCondition(Charge_conList));
   RCP<ConditionVisualDependency> Elec_pot1_Dep = rcp(
          new ConditionVisualDependency(Charge_andCon, ChargedFluid_List->getEntryRCP("CF6.2: Elec_pot_1"), true));

       /* end of set up for two dependees for the electrostatic condition */

      /*****************************************/
      /* add the dependencies for this section.*/
      /*****************************************/

   depSheet_Tramonto->addDependency(IcompVis_Dep); 
   depSheet_Tramonto->addDependency(IpolcompVis_Dep);
   depSheet_Tramonto->addDependency(DcoefIpolcompVis_Dep);
   depSheet_Tramonto->addDependency(DcoefIcompVis_Dep);
   depSheet_Tramonto->addDependency(Rhoicomp_Dep);
   depSheet_Tramonto->addDependency(Rhoipolcomp_Dep);
   depSheet_Tramonto->addDependency(GradDim_Dep);
   depSheet_Tramonto->addDependency(DensityUnit_Dep);
   depSheet_Tramonto->addDependency(Lconstrain_Dep);
   depSheet_Tramonto->addDependency(Diffusion_Dep);      
   depSheet_Tramonto->addDependency(DensityLength_Dep);
   depSheet_Tramonto->addDependency(DensityLength_DepPOL);
   depSheet_Tramonto->addDependency(ChargeParams_Dep);
   depSheet_Tramonto->addDependency(DielecPore_Dep);
   depSheet_Tramonto->addDependency(DielecConstUnit_Dep);
   depSheet_Tramonto->addDependency(Elec_pot1_Dep); 

   return;
}
/****************************************************************************************************************/
