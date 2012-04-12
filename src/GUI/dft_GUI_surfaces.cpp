using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_surfaces(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List)
{

     /* if we allow 10 surface geometry lists.  Need to construct all possible options for one list.  Then call this function multiple times for the different SGeomID_Lists.
        Here are the parameters we need and details of what we need to do:
        (1) Nwall_types, Nwall, Nlink, Lauto_center, Lauto_size --- easy --- surface list
        (2) Xtest_reflect ---- pretty easy  ---- surface list
        (3) Surface_type - validators are different for different Ndim (done)  --- surface geometry list.
        (4) Orientation - show only for appropriate surface_types ( and condition with Ndim )
        (5) WallParam - radius or halfwidth depending on surface_type (and condition with Ndim )
        (6) WallParam2 - half width or halflength depending on surface_type (and condition with Ndim)
        (7) WallParam3 - half width (Ndim=3 and cube surface)
        (7) Will there be any offsets applied to geometry parameters (random roughness, periodic or linear functions?)
        (8) Apply offsets to WallParam1, 2, or 3 (make this surface rough, periodic, or superimpose a linear function ?)
        (9) parameters for roughness */ 
        

     /* VALIDATORS*/

    RCP<EnhancedNumberValidator<int> > NwallValidator = rcp(new EnhancedNumberValidator<int>(0,1000,1));

   /* Set Surface PARAMETERS */
     Surface_List->set("S1: Number of Surfaces", 0, "Number of surfaces (or total surface subunits) in the problem",NwallValidator);
     Surface_List->set("S2: Number of macro surfaces",  Surface_List->get("S1: Number of Surfaces", (int)NULL), "Indicates if groups of surface form macrosurfaces.\n For example atoms that are part of a single molecule would have Nsurf_macro=1 and Nsurf>1)",NwallValidator);
     Surface_List->set("S3: Number of surface types", 1, "Number of different types surfaces (or total surface subunits) in the problem",NwallValidator);

  /* DEPENDENCIES */
     RCP<NumberVisualDependency<int> > NsurfType_Dep = rcp(
        new NumberVisualDependency<int>( Surface_List->getEntryRCP("S1: Number of Surfaces"), Surface_List->getEntryRCP("S3: Number of surface types")));

     RCP<NumberVisualDependency<int> > NsurfMacro_Dep = rcp(
        new NumberVisualDependency<int>( Surface_List->getEntryRCP("S1: Number of Surfaces"), Surface_List->getEntryRCP("S2: Number of macro surfaces")));

     /* DEPENDENCY SHEET ENTRIES*/
     depSheet_Tramonto->addDependency(NsurfType_Dep);
     depSheet_Tramonto->addDependency(NsurfMacro_Dep);

  return;
}

