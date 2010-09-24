using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_functionals(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List )
{
    /****************************************************************************************************************/
  /****************************** FUNCTIONAL CONTROL PARAMETER SECTION ********************************************/
  /****************************************************************************************************************/
    /**************************************/
    /* Define validators for this section.*/
    /**************************************/
    RCP<StringValidator> HSFuncValidator = rcp(
           new StringValidator(tuple<std::string>("Ideal Fluid / No Volume Exclusion","HS_Rosenfeld_Original", "HS_Rosenfeld_Zerocrossover","HS_White_Bear","HS_White_Bear2")));

    RCP<StringValidator> FuncPairPotValidator = rcp(
           new StringValidator(tuple<std::string>("No Mean Field Functional","UATT_core=[Umin; r=0,Rmin]","UATT_core=[0;r=0,Sigma]",
                "UATT_core=[U(sigma);r=0,Sigma]","UATT_core=[(0;r=0,Sigma),(Umin;r=Sigma,Rmin)]","UATTcore=[0;r=0,Rzero]")));

    RCP<StringValidator> CHARGEValidator = rcp(
           new StringValidator(tuple<std::string>("No Charge or No Poisson","Charge_Mean_Field","Charge_with_DeltaC_RPM", "Charge_with_DeltaC_General","Charge(MF)_with_Polarization")));

    RCP<StringValidator> POLYValidator = rcp(
           new StringValidator(tuple<std::string>("No Bonded Fluids","Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT","Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)","Polymer_JDC_iSAFT(comp)")));

    RCP<StringValidator> CalcTypeValidator = rcp(
           new StringValidator(tuple<std::string>("Equilibrium (homogeneous boundary conditions)","Equilibrium (inhomogeneous boudary conditions)","Steady State Diffusion (inhomogeneous boundaries)")));

    /***************************************************************/
    /* set up independent parameters that do not have dependencies */
    /***************************************************************/
    Functional_List->set("F0_Type_of_Calculation", "Equilibrium (homogeneous boundary conditions)", "Indicate the basic type of calculation to be performed.", CalcTypeValidator);
    Functional_List->set("F1_HS_Functional", "Ideal Fluid / No Volume Exclusion", "Type of Hard Sphere Functional (set to NONE for ideal gas or Poisson-Boltzman electrolyte)", HSFuncValidator);
    Functional_List->set("F2_PAIRPOTcore_Functional", "No Mean Field Functional", "Type of Core treatment for an extended pair potential \n using a strict mean field functional", FuncPairPotValidator);
    Functional_List->set("F3_CHARGE_Functional", "No Charge or No Poisson", "Functional for Charged systems (set to NONE if there are no charges in the system)", CHARGEValidator);
    Functional_List->set("F4_POLYMER_Functional", "No Bonded Fluids", "Functional for treating bonded systems", POLYValidator);

    /*******************************/
    /* define dependent parameters */
    /*******************************/

    /**********************************************************************************************/
    /* show the dependent parameters only if the independent parameters has a particular setting. */
    /**********************************************************************************************/

    /*****************************************/
    /* add the dependencies for this section.*/
    /*****************************************/

  /****************************************************************************************************************/
  /****************************** END FUNCTIONAL CONTROL PARAMETER SECTION ****************************************/
  /****************************************************************************************************************/
 
  return;
}

