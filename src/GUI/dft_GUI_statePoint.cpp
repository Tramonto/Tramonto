using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_StatePoint( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                  Teuchos::RCP<Teuchos::ParameterList> StatePoint_List,
                  Teuchos::RCP<Teuchos::ParameterList> Diffusion_List,
                  Teuchos::RCP<Teuchos::ParameterList> ChargedFluid_List)
{
  int idim;


        /**************************************/
        /* Define validators for this section.*/
        /**************************************/

    RCP<EnhancedNumberValidator<int> > DimValidator = rcp(new EnhancedNumberValidator<int>(0,2,1));
    RCP<StringValidator> DielecType_Validator = rcp(
           new StringValidator(tuple<std::string>("1 Dielec Const: uniform everywhere",
                       "2-Dielec Const: fluid and wall regions",
                       "3-Dielec Const: bulk fluid, fluid near wall, wall regions")));
 

        /*********************/
        /* set up parameters */
        /*********************/

   StatePoint_List->set("0: Temperature(K)", 300.0, "Enter the Temperature in Kelvin");
   StatePoint_List->set("1: Dimensionless Density Entry?", true, "Dimensionless densities are given in units of rho*sigma_ref^3.\n Set to false to enter in any other units (molar,g/cc etc).\n For dimensioned entry you must provide a conversion factor, Fac, where rho*sigma*3=(your densities)/Fac.");
   StatePoint_List->set("2: Density Conversion Factor",1.0,"Provide a conversion factor, Fac, where rho*sigma*3=(your densities)/Fac.");

   Array<double> RhoBulk0_Array( (Fluid_List->get<int>("F1_Ncomp")),0.0);
   StatePoint_List->set("3: Rho_b_0[icomp]", RhoBulk0_Array, "Set the bulk densities (homogeneous system), or the densities on the negative boundary (inhomogenous system)");

   Array<double> RhoBulk1_Array( (Fluid_List->get<int>("F1_Ncomp")),0.0);
   StatePoint_List->set("4: Rho_b_1[icomp]", RhoBulk1_Array, "Set the densities on the positive boundary  for an inhomogenous system");
   
   StatePoint_List->set("5: Direction of Gradient", 0, "Indicate direction where inhomogeneous boundaries should be applied (0=x,1=y,2=z)",DimValidator);
   StatePoint_List->set("6: Constrained Interface?", false, "Set to true to set rho[0]=0.5(rho[0]_left+rho[0]_right)\n at the midpoint of the domain. This approach can help to keep the interface from moving in a free interface calculation.");



   Array<double> DiffCoeff_Array( (Fluid_List->get<int>("F1_Ncomp")),0.0);
   Diffusion_List->set("D1: Diff_Coeff[icomp]", DiffCoeff_Array, "Set the diffusion coefficientf sof each component.");

   Diffusion_List->set("D2: distance for constant Mu region", 0.0, "Set the distance from the computational boundaries\n in the direction of the chemical potential gradient\n where the chemical potential is to be held constant.\n Facilitates comparison with GCMD calculations.");

   Diffusion_List->set("D3: bulk velocity", 0.0, "Set a bulk velocity term in diffusion calculation.");


   ChargedFluid_List->set("CF1: Type Dielectric Constant(s)","1 Dielec Const: uniform everywhere","Indicate type of dielectric constants to be entered for the problem.",DielecType_Validator);
   ChargedFluid_List->set("CF2: Entry of Relative Dielectric Constant(s)?",true,"true indicates that dielectric constants will be entered as D/D_reference\n (for example, D_ref may be the dielectric constant of the bulk solvent).\n If set to false, you will be asked for provide a reference value, D_reference for computation of relative dielectric constants.");
   ChargedFluid_List->set("CF3: Reference Dielectric Constant",78.5,"Enter a reference dielectric constant for calculation of relative values D/D_reference.\n Note that 78.5 is a typical value used for water");

   ChargedFluid_List->set("CF4: Dielec Const Bulk Fluid",1.0,"Enter the dielectric constant of the bulk fluid.");
   ChargedFluid_List->set("CF5: Dielec Const Near Wall Fluid",1.0,"Enter the dielectric constant of the fluid in the pore or near the wall.");
   ChargedFluid_List->set("CF6: Size of Near Wall region",1.0,"Enter the distance from the wall to be considered in the Near Wall region.");

   ChargedFluid_List->set("CF7: Plasma Parameter",1.0,"Enter the dimensionless plasma parameter defined as:\n\t e^2/(4pi kT[D*epsilon_0*sigma_ref]\n where k=1.3807e-23J/K; epsilon_0=8.85419e-12C^2/(Jm); e=1.60219e-19C; pi=3.14159.\n T is the temperature, D is a reference dielectric constant, and sigma is a reference length.");

   

   ChargedFluid_List->set("CF8.0: Dimensionless Entry of Electrostatic Potential(s)?",true,"true indicates that electrostatic potentials will be entered as e(phi)/kT.\n If set to false, you will be required to enter a temperature in Kelvin units.");

   ChargedFluid_List->set("CF8.1: Elec_pot_0", 0.0, "Set a bulk electrostatic potential (homogeneous boundary case), or \n set the bulk electrostatic potential on the negative boundary (inhomogeneous boundary or diffusion).");

   ChargedFluid_List->set("CF8.2: Elec_pot_1", 0.0, "Set the bulk electrostatic potential on the positive boundary (inhomogeneous boundary or diffusion).");


        /************************/
        /* set up dependencies */
        /************************/
   RCP<StringVisualDependency> GradDim_Dep = rcp(
       new StringVisualDependency( Functional_List->getEntryRCP("F0_Type_of_Calculation"),
                                   StatePoint_List->getEntryRCP("5: Direction of Gradient"), 
                                   tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)",
                                   "Steady State Diffusion (inhomogeneous boundaries)")));

       /* set up two dependees for the temperature condition ...note there used to be a false here so that a negative entry meant enter dimensioned units...may need to work on this*/
   RCP<BoolCondition> TempCon1=rcp(new BoolCondition(Fluid_List->getEntryRCP("F5_Dimensionless_Energy_Entry")));
   RCP<BoolCondition> TempCon2=rcp(new BoolCondition(ChargedFluid_List->getEntryRCP("CF8.0: Dimensionless Entry of Electrostatic Potential(s)?")));
   Condition::ConstConditionList Temp_conList = tuple<RCP<const Condition> >(TempCon1,TempCon2);
   RCP<OrCondition>Temp_orCon=rcp(new OrCondition(Temp_conList));
   RCP<ConditionVisualDependency> EnergyUnit_Dep = rcp(new ConditionVisualDependency(Temp_orCon, StatePoint_List->getEntryRCP("0: Temperature(K)")));
       /* end of set up two dependees for the temperature condition */

   RCP<BoolVisualDependency> DensityUnit_Dep = rcp( 
        new BoolVisualDependency( StatePoint_List->getEntryRCP("1: Dimensionless Density Entry?"),
                                  StatePoint_List->getEntryRCP("2: Density Conversion Factor"),false));

   RCP<StringVisualDependency> Lconstrain_Dep = rcp(
       new StringVisualDependency( Functional_List->getEntryRCP("F0_Type_of_Calculation"),
                                   StatePoint_List->getEntryRCP("6: Constrained Interface?"),
                                   tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)")));

   Dependency::ParameterEntryList DensityLengthDependents;
   DensityLengthDependents.insert(StatePoint_List->getEntryRCP("3: Rho_b_0[icomp]"));
   DensityLengthDependents.insert(StatePoint_List->getEntryRCP("4: Rho_b_1[icomp]"));
   RCP<NumberArrayLengthDependency<int,double> > DensityLength_Dep = rcp(
           new NumberArrayLengthDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"), DensityLengthDependents));

   RCP<StringVisualDependency> RhoB1_Dep = rcp(
       new StringVisualDependency( Functional_List->getEntryRCP("F0_Type_of_Calculation"),
                                   StatePoint_List->getEntryRCP("4: Rho_b_1[icomp]"), 
                                   tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)",
                                   "Steady State Diffusion (inhomogeneous boundaries)")));

   Dependency::ParameterEntryList DiffusionDependents;
   DiffusionDependents.insert(Diffusion_List->getEntryRCP("D1: Diff_Coeff[icomp]"));
   DiffusionDependents.insert(Diffusion_List->getEntryRCP("D2: distance for constant Mu region"));
   DiffusionDependents.insert(Diffusion_List->getEntryRCP("D3: bulk velocity"));

   RCP<StringVisualDependency> Diffusion_Dep = rcp(
       new StringVisualDependency( Functional_List->getEntryRCP("F0_Type_of_Calculation"),DiffusionDependents,
           tuple<std::string>("Steady State Diffusion (inhomogeneous boundaries)")));

   Dependency::ParameterEntryList ChargeParamsDependents;
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF1: Type Dielectric Constant(s)"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF2: Entry of Relative Dielectric Constant(s)?"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF4: Dielec Const Bulk Fluid"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF7: Plasma Parameter"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF8.0: Dimensionless Entry of Electrostatic Potential(s)?"));
   ChargeParamsDependents.insert(ChargedFluid_List->getEntryRCP("CF8.1: Elec_pot_0"));

   RCP<StringVisualDependency> ChargeParams_Dep = rcp(
       new StringVisualDependency( Functional_List->getEntryRCP("F3_CHARGE_Functional"),ChargeParamsDependents,
           tuple<std::string>("Charge_Mean_Field","Charge_with_DeltaC_RPM","Charge_with_DeltaC_General","Charge(MF)_with_Polarization")));

   RCP<BoolVisualDependency> DielecConstUnit_Dep = rcp( 
        new BoolVisualDependency( ChargedFluid_List->getEntryRCP("CF2: Entry of Relative Dielectric Constant(s)?"),
                                  ChargedFluid_List->getEntryRCP("CF3: Reference Dielectric Constant"),false));

   Dependency::ParameterEntryList DielecPoreDependents;
   DielecPoreDependents.insert(ChargedFluid_List->getEntryRCP("CF5: Dielec Const Near Wall Fluid"));
   DielecPoreDependents.insert(ChargedFluid_List->getEntryRCP("CF6: Size of Near Wall region"));
   RCP<StringVisualDependency> DielecPore_Dep = rcp(
       new StringVisualDependency( ChargedFluid_List->getEntryRCP("CF1: Type Dielectric Constant(s)"),DielecPoreDependents,
           tuple<std::string>("2-Dielec Const: fluid and wall regions","3-Dielec Const: bulk fluid, fluid near wall, wall regions")));

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
          new ConditionVisualDependency(Charge_andCon, ChargedFluid_List->getEntryRCP("CF8.2: Elec_pot_1"), true));

       /* end of set up for two dependees for the electrostatic condition */


      /*****************************************/
      /* add the dependencies for this section.*/
      /*****************************************/

   depSheet_Tramonto->addDependency(GradDim_Dep);
   depSheet_Tramonto->addDependency(EnergyUnit_Dep);
   depSheet_Tramonto->addDependency(DensityUnit_Dep);
   depSheet_Tramonto->addDependency(Lconstrain_Dep);
   depSheet_Tramonto->addDependency(DensityLength_Dep);
   depSheet_Tramonto->addDependency(RhoB1_Dep);
   depSheet_Tramonto->addDependency(Diffusion_Dep);
   depSheet_Tramonto->addDependency(ChargeParams_Dep);
   depSheet_Tramonto->addDependency(DielecPore_Dep);
   depSheet_Tramonto->addDependency(DielecConstUnit_Dep);
   depSheet_Tramonto->addDependency(Elec_pot1_Dep);


  return;
}

