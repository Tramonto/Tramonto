/*#include "Teuchos_ParameterList.hpp"*/
using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_mesh( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Mesh_List) 
{
  int idim;

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
    Mesh_List->set("M1_Ndim", 1, "Number of dimensions for numerical computation.", NdimValidator);

    Mesh_List->set("M2_Dimensionless_Distance_Entry", true, "True if all length variables will be entered in dimensionless units (L/sigma).\n Otherwise to enter distance parameters in some other units (Angstroms, nm etc),\n you will be asked to provide a conversion factor such that L/sigma=(your length units)/Fac.");

    /* define dependent parameters */
    Mesh_List->set("M3_Reference_Length", 1.0, "Provide the unit conversion factor such that L/sigma=(your length units)/Reference_Length");

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
           new BoolVisualDependency( 
               Mesh_List->getEntryRCP("M2_Dimensionless_Distance_Entry"),
               Mesh_List->getEntryRCP("M3_Reference_Length"),
               false)
    );

    RCP<NumberArrayLengthDependency<int, double> > EsizeArrayDep = rcp(
          new NumberArrayLengthDependency<int, double> ( 
              Mesh_List->getEntryRCP("M1_Ndim"),
              Mesh_List->getEntryRCP("M5_Esize_x"))
    );

    RCP<NumberArrayLengthDependency<int, double> > SizeArrayDep = rcp(
          new NumberArrayLengthDependency<int, double> ( 
              Mesh_List->getEntryRCP("M1_Ndim"),
              Mesh_List->getEntryRCP("M4_Size_x"))
    );

    RCP<NumberArrayLengthDependency<int, std::string> > BCnegArrayDep = rcp(
          new NumberArrayLengthDependency<int, std::string> ( 
              Mesh_List->getEntryRCP("M1_Ndim"), 
              Mesh_List->getEntryRCP("M6_Boundary_Conditions_XMIN"))
    );

    RCP<NumberArrayLengthDependency<int, std::string> > BCposArrayDep = rcp(
          new NumberArrayLengthDependency<int, std::string> ( 
              Mesh_List->getEntryRCP("M1_Ndim"), 
              Mesh_List->getEntryRCP("M7_Boundary_Conditions_XMAX"))
    );


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

