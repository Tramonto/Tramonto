/*#include "Teuchos_ParameterList.hpp"*/
using namespace std;
#include <iostream>
/***** push these to common dft_GUI.h file *****
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"	
#include "Teuchos_Version.hpp"
#include "Optika_GUI.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"
**************************************************/
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

        /*********************/
        /* set up parameters */
        /*********************/

   StatePoint_List->set("SP1: Direction of Gradient", 0, "Indicate direction where inhomogeneous boundaries should be applied (0=x,1=y,2=z)",DimValidator);
   StatePoint_List->set("SP2: Constrained Interface?", false, "Set to true to set rho[0]=0.5(rho[0]_left+rho[0]_right)\n at the midpoint of the domain. This approach can help to keep the interface from moving in a free interface calculation.");

   Array<double> RhoBulk0_Array( (Fluid_List->get<int>("F1_Ncomp")),0.0);
   StatePoint_List->set("SP3: Rho_b_0[icomp]", RhoBulk0_Array, "Set the bulk densities (homogeneous system), or the densities on the negative boundary (inhomogenous system)");

   Array<double> RhoBulk1_Array( (Fluid_List->get<int>("F1_Ncomp")),0.0);
   StatePoint_List->set("SP4: Rho_b_1[icomp]", RhoBulk1_Array, "Set the densities on the positive boundary  for an inhomogenous system");

   Array<double> DiffCoeff_Array( (Fluid_List->get<int>("F1_Ncomp")),0.0);
   Diffusion_List->set("D1: Diff_Coeff[icomp]", DiffCoeff_Array, "Set the diffusion coefficientf sof each component.");

   Diffusion_List->set("D2: distance for constant Mu region", 0.0, "Set the distance from the computational boundaries\n in the direction of the chemical potential gradient\n where the chemical potential is to be held constant.\n Facilitates comparison with GCMD calculations.");

   Diffusion_List->set("D3: bulk velocity", 0.0, "Set a bulk velocity term in diffusion calculation.");

   ChargedFluid_List->set("CF1: Elec_pot_0", 0.0, "Set a bulk electrostatic potential (homogeneous boundary case), or \n set the bulk electrostatic potential on the negative boundary (inhomogeneous boundary or diffusion).");

   ChargedFluid_List->set("CF2: Elec_pot_1", 0.0, "Set the bulk electrostatic potential on the positive boundary (inhomogeneous boundary or diffusion).");


        /************************/
        /* set up dependencies */
        /************************/
   RCP<StringVisualDependency> GradDim_Dep = rcp(
       new StringVisualDependency( "F0_Type_of_Calculation",Functional_List,"SP1: Direction of Gradient", StatePoint_List,
           tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)","Steady State Diffusion (inhomogeneous boundaries)")));

   RCP<StringVisualDependency> Lconstrain_Dep = rcp(
       new StringVisualDependency( "F0_Type_of_Calculation",Functional_List,"SP2: Constrained Interface?", StatePoint_List,
           tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)")));

   RCP<NumberArrayLengthDependency> RhoB0Length_Dep = rcp(
           new NumberArrayLengthDependency( "F1_Ncomp", Fluid_List, "SP3: Rho_b_0[icomp]",StatePoint_List));

   RCP<NumberArrayLengthDependency> RhoB1Length_Dep = rcp(
           new NumberArrayLengthDependency( "F1_Ncomp", Fluid_List, "SP4: Rho_b_1[icomp]",StatePoint_List));

   RCP<StringVisualDependency> RhoB1_Dep = rcp(
       new StringVisualDependency( "F0_Type_of_Calculation",Functional_List,"SP4: Rho_b_1[icomp]", StatePoint_List,
           tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)","Steady State Diffusion (inhomogeneous boundaries)")));

   RCP<StringVisualDependency> DiffCoeff_Dep = rcp(
       new StringVisualDependency( "F0_Type_of_Calculation",Functional_List,"D1: Diff_Coeff[icomp]", Diffusion_List,
           tuple<std::string>("Steady State Diffusion (inhomogeneous boundaries)")));

   RCP<StringVisualDependency> XconstMu_Dep = rcp(
       new StringVisualDependency( "F0_Type_of_Calculation",Functional_List,"D2: distance for constant Mu region", Diffusion_List,
           tuple<std::string>("Steady State Diffusion (inhomogeneous boundaries)")));

   RCP<StringVisualDependency> Velocity_Dep = rcp(
       new StringVisualDependency( "F0_Type_of_Calculation",Functional_List,"D3: bulk velocity", Diffusion_List,
           tuple<std::string>("Steady State Diffusion (inhomogeneous boundaries)")));

   RCP<StringVisualDependency> Elec_pot0_Dep = rcp(
       new StringVisualDependency( "F3_CHARGE_Functional",Functional_List,"CF1: Elec_pot_0", ChargedFluid_List,
           tuple<std::string>("Charge_Mean_Field","Charge_with_DeltaC_RPM","Charge_with_DeltaC_General","Charge(MF)_with_Polarization")));

       /* set up two dependees for the electrostatic condition */
   RCP<StringCondition> ChargeCon1 = rcp(new StringCondition("F3_CHARGE_Functional", Functional_List, 
           tuple<std::string>("Charge_Mean_Field","Charge_with_DeltaC_RPM","Charge_with_DeltaC_General","Charge(MF)_with_Polarization"),true));
 
   RCP<StringCondition> ChargeCon2 = rcp(new StringCondition("F0_Type_of_Calculation", Functional_List, 
           tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)","Steady State Diffusion (inhomogeneous boundaries)"),true));

   Condition::ConditionList Charge_conList = tuple<RCP<Condition> >(ChargeCon1, ChargeCon2);
   RCP<AndCondition> Charge_andCon = rcp(new AndCondition(Charge_conList));

   RCP<ConditionVisualDependency> Elec_pot1_Dep = rcp(new ConditionVisualDependency(Charge_andCon, "CF2: Elec_pot_1", ChargedFluid_List, true));
       /* end of set up for two dependees for the electrostatic condition */


      /*****************************************/
      /* add the dependencies for this section.*/
      /*****************************************/

   depSheet_Tramonto->addDependency(GradDim_Dep);
   depSheet_Tramonto->addDependency(Lconstrain_Dep);
   depSheet_Tramonto->addDependency(RhoB0Length_Dep);
   depSheet_Tramonto->addDependency(RhoB1Length_Dep);
   depSheet_Tramonto->addDependency(RhoB1_Dep);
   depSheet_Tramonto->addDependency(DiffCoeff_Dep);
   depSheet_Tramonto->addDependency(XconstMu_Dep);
   depSheet_Tramonto->addDependency(Velocity_Dep);
   depSheet_Tramonto->addDependency(Elec_pot0_Dep);
   depSheet_Tramonto->addDependency(Elec_pot1_Dep);


  return;
}

