using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_Diffusion( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                  Teuchos::RCP<Teuchos::ParameterList> StatePoint_List)
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


        /************************/
        /* set up dependencies */
        /************************/

   Dependency::ParameterParentMap InhomogeneousBCDependents;
   InhomogeneousBCDependents.insert(std::pair<std::string, RCP<ParameterList> >("SP1: Direction of Gradient", StatePoint_List));
   InhomogeneousBC.insert(std::pair<std::string, RCP<ParameterList> >("SP4: Rho_b_1[icomp]", StatePoint_List));

   RCP<StringVisualDependency> InhomogeneousBC_Dep = rcp(
       new StringVisualDependency( "F0_Type_of_Calculation",Functional_List,InhomogeneousBCDependents,
           tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)","Steady State Diffusion (inhomogeneous boundaries)")));

   RCP<StringVisualDependency> Lconstrain_Dep = rcp(
       new StringVisualDependency( "F0_Type_of_Calculation",Functional_List,"SP2: Constrained Interface?", StatePoint_List,
           tuple<std::string>("Equilibrium (inhomogeneous boudary conditions)")));

   RCP<NumberArrayLengthDependency<int,double> > RhoB0Length_Dep = rcp(
           new NumberArrayLengthDependency<int,double>( Fluid_List->getEntryRCP("F1_Ncomp"), StatePoint_List->getEntryRCP("SP3: Rho_b_0[icomp]")));

   RCP<NumberArrayLengthDependency<int,double> > RhoB1Length_Dep = rcp(
           new NumberArrayLengthDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"),StatePoint_List->getEntryRCP("SP4: Rho_b_1[icomp]")));

      /*****************************************/
      /* add the dependencies for this section.*/
      /*****************************************/

   depSheet_Tramonto->addDependency(InhomogeneousBC_Dep);
   depSheet_Tramonto->addDependency(Lconstrain_Dep);
   depSheet_Tramonto->addDependency(RhoB0Length_Dep);
   depSheet_Tramonto->addDependency(RhoB1Length_Dep);


  return;
}

