using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_surfaces_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List)
{
  /* VALIDATORS*/
  RCP<EnhancedNumberValidator<int> > NwallValidator = rcp(new EnhancedNumberValidator<int>(0,1000,1));
  RCP<EnhancedNumberValidator<int> > NwallTypeValidator = rcp(new EnhancedNumberValidator<int>(0,12,1));

  /* Set Surface PARAMETERS */
  Surface_List->set("S1: Number of Surfaces", 0, "Number of surfaces (or total surface subunits) in the problem",NwallValidator);
  Surface_List->set("S2: Number of macro surfaces",  Surface_List->get<int>("S1: Number of Surfaces"), "Indicates if groups of surface form macrosurfaces.\n For example atoms that are part of a single molecule would have Nsurf_macro=1 and Nsurf>1)",NwallValidator);
  Surface_List->set("S3: Number of surface types", 1, "Number of different types surfaces (or total surface subunits) in the problem",NwallTypeValidator);

   return;
}
/*******************************************************************************************/
void dft_GUI_surfaces_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List)
{
  /* VALIDATORS*/
  RCP<EnhancedNumberValidator<int> > NwallValidator = rcp(new EnhancedNumberValidator<int>(0,1000,1));
  RCP<EnhancedNumberValidator<int> > NwallTypeValidator = rcp(new EnhancedNumberValidator<int>(0,12,1));

  Surface_List->set("S1: Number of Surfaces", Nwall, "Number of surfaces (or total surface subunits) in the problem",NwallValidator);
  Surface_List->set("S2: Number of macro surfaces",  Nlink, "Indicates if groups of surface form macrosurfaces.\n For example atoms that are part of a single molecule would have Nsurf_macro=1 and Nsurf>1)",NwallValidator);
  Surface_List->set("S3: Number of surface types", Nwall_type, "Number of different types surfaces (or total surface subunits) in the problem",NwallTypeValidator);
  return;
}
/*******************************************************************************************/
void dft_GUI_surfaces_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List)
{
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
/*******************************************************************************************/

