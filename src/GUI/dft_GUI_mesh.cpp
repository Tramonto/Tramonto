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

void dft_GUI_mesh( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Mesh_List) 
{
  int idim;
  /****************************** DIMENSION PARAMETER SECTION **********************************************************/
  /* begin development of GUI using old input structure found in dft_input.dat.  The first part of the input
     file identifies what type of input we should expect - dimensionless or dimensioned parameters for several
     different variables - length, density, temperature, and dielectric constant. 
     To accomplish this, both validators and dependency lists are needed.  Also, all these parameters will be put in
     a sublist called Basic Parameter Definitions */



    /* define the sublist */

    RCP<ParameterList> Units_List = sublist(Tramonto_List,"Sect. 0: Units for Parameter Entry");

  /* Define validators for this section.*/

   /* set up independent parameters that do not have dependencies */

  Units_List->set("Dimensionless_Entry_Energy", true, "Will all energy variables be entered in dimensionless units (E/kT)?");
  Units_List->set("Dimensionless_Entry_Density", true, "Will all density variables be entered in dimensionless units (rho*sigma^3)?");
  /* set up dependent (hidden) parameters */
  Units_List->set("Reference_Temperature", 1.0, "Provide the temperature in Kelvin. \nTramonto will set up dimensionless energy variables using eps/(kT)");
  Units_List->set("Reference_Density", 1.0, "Provide a reference density in units consistent with your data entry. \n Reference_Density must be set such that the input density variables (rho)\n will be reduced to rho*sigma_r^3  when the input values are divided by Reference_density (rho sigma^3=rho/Reference_Density)");



        /* define dependent parameters */

        /* in this case the default value may be set to two different parameters - or in one case the user must set this parameter.*/

/*         cout<<"output from the Tramonto_List-get call is"<<Tramonto_List->get("ParamEntryType", (std::string)"")<<endl;*/

        /* show the dependent parameters only if the independent parameters has a particular setting. */


      RCP<BoolVisualDependency> DimensionDep_temp = rcp(
           new BoolVisualDependency( "Dimensionless_Entry_Energy",Units_List,"Reference_Temperature", Units_List, false));

      RCP<BoolVisualDependency> DimensionDep_density = rcp(
           new BoolVisualDependency( "Dimensionless_Entry_Density",Units_List,"Reference_Density", Units_List, false));

      /*
       * add the dependencies for this section.
       */
       depSheet_Tramonto->addDependency(DimensionDep_density);
       depSheet_Tramonto->addDependency(DimensionDep_temp);


  /****************************************************************************************************************/
  /****************************** END DIMENSION PARAMETER SECTION *************************************************/
  /****************************************************************************************************************/

  /****************************************************************************************************************/
  /****************************** MESH PARAMETER SECTION **********************************************************/
  /****************************************************************************************************************/
  /* This section creates parameters for: the size of the computational domain, the mesh density, and the boundary conditions */


    /* define the sublist */
/*    RCP<ParameterList> Mesh_List = sublist(Tramonto_List,"Sect. 1: Computational Domain Parameters");*/

    /* Define validators for this section.*/
    RCP<EnhancedNumberValidator<int> > NdimValidator = rcp(new EnhancedNumberValidator<int>(1,3,1));

    RCP<StringValidator> BCValidator = rcp(
           new StringValidator(tuple<std::string>("Bulk","In_Wall","Reflect","Periodic","Last_node","Last_node_restart")));
    RCP<ArrayStringValidator> BCArrayValidator = rcp(new ArrayStringValidator(BCValidator));

    /* set up independent parameters that do not have dependencies */
    Mesh_List->set("M2_Dimensionless_Entry_Mesh_TF", true, "Will all length variables be entered in dimensionless units (L/sigma)?");
    Mesh_List->set("M1_Ndim", 1, "Number of dimensions for numerical computation.", NdimValidator);

    /* define dependent parameters */
    Mesh_List->set("M3_Reference_Length", 1.0, "Provide a reference length (sigma) in the same units \n as you wish to use for data entry.  Tramonto will set up \n dimensionless parameters as param/Reference_length.");

    Array<double> Size_x_Array( Mesh_List->get<int>("M1_Ndim"),10.0);
    Mesh_List->set("M4_Size_x", Size_x_Array, "the size of the computational domain : Ndim entries");

        /*cout << "Ndim="<<Mesh_List->get<int>("M1_Ndim")<<endl;

    for (idim=0; idim<Mesh_List->get<int>("M1_Ndim"); idim++){
        cout << "idim="<<idim<<"Size_x="<<Mesh_List->get<double>("M4_Size_x")<<endl;
    }*/

    Array<double> Esize_x_Array( Mesh_List->get<int>("M1_Ndim"),0.1);
    Mesh_List->set("M5_Esize_x", Esize_x_Array, "the mesh spacing for a cartesian mesh : Ndim entries");

    Array<std::string> BCneg_Array( Mesh_List->get<int>("M1_Ndim"),"Bulk");
    Mesh_List->set("M6_Boundary_Conditions_XMIN", BCneg_Array, "boundary condition to apply at X_min,Y_min,Z_min: Ndim entries",BCArrayValidator);

    Array<std::string> BCpos_Array( Mesh_List->get<int>("M1_Ndim"),"Bulk");
    Mesh_List->set("M7_Boundary_Conditions_XMAX", BCpos_Array, "boundary condition to apply for X_max,Y_max,Z_mix: Ndim entries",BCArrayValidator);


    /* show the dependent parameters only if the independent parameters has a particular setting. */
     RCP<BoolVisualDependency> UnitEntryDep_mesh = rcp(
           new BoolVisualDependency( "M2_Dimensionless_Entry_Mesh_TF",Mesh_List,"M3_Reference_Length", Mesh_List, false));

    RCP<NumberArrayLengthDependency> EsizeArrayDep = rcp(
          new NumberArrayLengthDependency( "M1_Ndim", Mesh_List, "M5_Esize_x",Mesh_List));

    RCP<NumberArrayLengthDependency> SizeArrayDep = rcp(
          new NumberArrayLengthDependency( "M1_Ndim", Mesh_List, "M4_Size_x",Mesh_List));

    RCP<NumberArrayLengthDependency> BCnegArrayDep = rcp(
          new NumberArrayLengthDependency( "M1_Ndim", Mesh_List, "M6_Boundary_Conditions_XMIN",Mesh_List));

    RCP<NumberArrayLengthDependency> BCposArrayDep = rcp(
          new NumberArrayLengthDependency( "M1_Ndim", Mesh_List, "M7_Boundary_Conditions_XMAX",Mesh_List));


    /* add the dependencies for this section.*/
       depSheet_Tramonto->addDependency(UnitEntryDep_mesh);
       depSheet_Tramonto->addDependency(SizeArrayDep);
       depSheet_Tramonto->addDependency(EsizeArrayDep);
       depSheet_Tramonto->addDependency(BCnegArrayDep);
       depSheet_Tramonto->addDependency(BCposArrayDep);

  /****************************************************************************************************************/
  /****************************** END MESH PARAMETER SECTION ******************************************************/
  /****************************************************************************************************************/
   
  /****************************************************************************************************************/
  /****************************** ANY NEW SUBSECTION PARAMETER SECTION ********************************************/
  /****************************************************************************************************************/
    /* SUBLIST(S) */
    /* DEFINE VALIDATORS*/ 
    /* INDEPENDENT PARAMETERS */
    /* DEPENDENT PARAMETERS */
    /* DEPENDENCIES */
    /* DEPENDENCY SHEET ENTRIES*/
  /****************************************************************************************************************/
  /****************************** END ANY NEW SUBSECTION PARAMETER SECTION ****************************************/
  /****************************************************************************************************************/

  return;
}

