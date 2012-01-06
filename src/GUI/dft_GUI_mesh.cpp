/*#include "Teuchos_ParameterList.hpp"*/
using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_mesh( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Mesh_List) 
{
  int idim;
  bool set_defaults_from_old_format_file=true;
  string BC_array_tmp[3];

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

    /* set up all parameters for this section */




    if (set_defaults_from_old_format_file){
        Mesh_List->set("M1_Ndim", Ndim, "Number of dimensions for numerical computation.", NdimValidator);

        if (Length_ref<=0.) Mesh_List->set("M2_Dimensionless_Distance_Entry", true, "True if all length variables will be entered in dimensionless units (L/sigma).\n Otherwise to enter distance parameters in some other units (Angstroms, nm etc),\n you will be asked to provide a conversion factor such that L/sigma=(your length units)/Fac.");
         else Mesh_List->set("M2_Dimensionless_Distance_Entry", false, "True if all length variables will be entered in dimensionless units (L/sigma).\n Otherwise to enter distance parameters in some other units (Angstroms, nm etc),\n you will be asked to provide a conversion factor such that L/sigma=(your length units)/Fac.");

         Mesh_List->set("M3_Reference_Length", Length_ref, "Provide the unit conversion factor such that L/sigma=(your length units)/Reference_Length");

         Array<double> Size_x_Array(Size_x, Size_x+Ndim);
         Mesh_List->set("M4_Size_x", Size_x_Array, "the size of the computational domain : Ndim entries");

         Array<double> Esize_x_Array(Esize_x, Esize_x+Ndim);
         Mesh_List->set("M5_Esize_x", Esize_x_Array, "the mesh spacing for a cartesian mesh : Ndim entries");

         for (idim=0;idim<Ndim;idim++){
           if (Type_bc[idim][0]==IN_BULK) BC_array_tmp[idim]="Bulk";
           else if (Type_bc[idim][0]==IN_WALL) BC_array_tmp[idim]="In_Wall" ;
           else if (Type_bc[idim][0]==REFLECT)  BC_array_tmp[idim]="Reflect";
           else if (Type_bc[idim][0]==PERIODIC)  BC_array_tmp[idim]="Periodic";
           else if (Type_bc[idim][0]==LAST_NODE) BC_array_tmp[idim]="Last_node";
           else if (Type_bc[idim][0]==LAST_NODE_RESTART) BC_array_tmp[idim]="Last_node_restart";
         }
         Array<std::string> BCneg_Array(BC_array_tmp,BC_array_tmp+Ndim);
         Mesh_List->set("M6_Boundary_Conditions_XMIN", BCneg_Array, "boundary condition to apply at X_min,Y_min,Z_min: Ndim entries",BCArrayValidator);

         for (idim=0;idim<Ndim;idim++){
           if (Type_bc[idim][1]==IN_BULK) BC_array_tmp[idim]="Bulk";
           else if (Type_bc[idim][1]==IN_WALL) BC_array_tmp[idim]="In_Wall" ;
           else if (Type_bc[idim][1]==REFLECT)  BC_array_tmp[idim]="Reflect";
           else if (Type_bc[idim][1]==PERIODIC)  BC_array_tmp[idim]="Periodic";
           else if (Type_bc[idim][1]==LAST_NODE) BC_array_tmp[idim]="Last_node";
           else if (Type_bc[idim][1]==LAST_NODE_RESTART) BC_array_tmp[idim]="Last_node_restart";
         }
         Array<std::string> BCpos_Array(BC_array_tmp,BC_array_tmp+Ndim);
         Mesh_List->set("M7_Boundary_Conditions_XMAX", BCpos_Array, "boundary condition to apply for X_max,Y_max,Z_mix: Ndim entries",BCArrayValidator);
    }
    else{
        Mesh_List->set("M1_Ndim", 1, "Number of dimensions for numerical computation.", NdimValidator);
    
        Mesh_List->set("M2_Dimensionless_Distance_Entry", true, "True if all length variables will be entered in dimensionless units (L/sigma).\n Otherwise to enter distance parameters in some other units (Angstroms, nm etc),\n you will be asked to provide a conversion factor such that L/sigma=(your length units)/Fac.");

        Mesh_List->set("M3_Reference_Length", 1.0, "Provide the unit conversion factor such that L/sigma=(your length units)/Reference_Length");

        Array<double> Size_x_Array( Mesh_List->get<int>("M1_Ndim"),10.0);
        Mesh_List->set("M4_Size_x", Size_x_Array, "the size of the computational domain : Ndim entries");

        Array<double> Esize_x_Array( Mesh_List->get<int>("M1_Ndim"),0.1);
        Mesh_List->set("M5_Esize_x", Esize_x_Array, "the mesh spacing for a cartesian mesh : Ndim entries");

       Array<std::string> BCneg_Array( Mesh_List->get<int>("M1_Ndim"),"Bulk");
       Mesh_List->set("M6_Boundary_Conditions_XMIN", BCneg_Array, "boundary condition to apply at X_min,Y_min,Z_min: Ndim entries",BCArrayValidator);

       Array<std::string> BCpos_Array( Mesh_List->get<int>("M1_Ndim"),"Bulk");
       Mesh_List->set("M7_Boundary_Conditions_XMAX", BCpos_Array, "boundary condition to apply for X_max,Y_max,Z_mix: Ndim entries",BCArrayValidator);
    }

    
    /* show the dependent parameters only if the independent parameters has a particular setting. */
    RCP<BoolVisualDependency> UnitEntryDep_mesh = rcp(
           new BoolVisualDependency(Mesh_List->getEntryRCP("M2_Dimensionless_Distance_Entry"),Mesh_List->getEntryRCP("M3_Reference_Length"),false));

    RCP<NumberArrayLengthDependency<int,double> > EsizeArrayDep = rcp(
          new NumberArrayLengthDependency<int,double>(Mesh_List->getEntryRCP("M1_Ndim"),Mesh_List->getEntryRCP("M5_Esize_x")));

    RCP<NumberArrayLengthDependency<int,double> > SizeArrayDep = rcp(
          new NumberArrayLengthDependency<int,double>( Mesh_List->getEntryRCP("M1_Ndim"), Mesh_List->getEntryRCP("M4_Size_x")));

    RCP<NumberArrayLengthDependency<int,string> > BCnegArrayDep = rcp(
          new NumberArrayLengthDependency<int,string>( Mesh_List->getEntryRCP("M1_Ndim"), Mesh_List->getEntryRCP("M6_Boundary_Conditions_XMIN")));

    RCP<NumberArrayLengthDependency<int,string> > BCposArrayDep = rcp(
          new NumberArrayLengthDependency<int,string>( Mesh_List->getEntryRCP("M1_Ndim"), Mesh_List->getEntryRCP("M7_Boundary_Conditions_XMAX")));


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

