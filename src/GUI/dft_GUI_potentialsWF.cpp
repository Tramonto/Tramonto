using namespace std;
#include <iostream>
/*#include "dft_GUI.h"*/
#include "dft_globals_const.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

/***********************************************************************************************************************************/
void dft_GUI_potentialsWF_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Surface_List, 
                         Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsWF_List )
{
    /********************/
    /* setup parameters */
    /********************/

    TwoDArray<double> SigmaWF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Surface_List->get<int>("S3: Number of surface types")),1.0);
    PotentialsWF_List->set("WF1: Sigma_wf", SigmaWF_Array, "Sigma - Characteristic diameter fluid-fluid pair interactions");

    TwoDArray<double> EpsWF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Surface_List->get<int>("S3: Number of surface types")),1.0);
    PotentialsWF_List->set("WF2: Eps_wf", EpsWF_Array, "Eps - Energy prefactor for fluid-fluid pair interactions");

    TwoDArray<double> CutWF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Surface_List->get<int>("S3: Number of surface types")),3.0);
    PotentialsWF_List->set("WF3: Cut_wf", CutWF_Array, "rcut - Cutoff distance for fluid-fluid pair interactions");

    TwoDArray<double> EpsYukawaWF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Surface_List->get<int>("S3: Number of surface types")),1.0);
    PotentialsWF_List->set("WF4: EpsYukawa_wf", EpsYukawaWF_Array, "AYuk - Energy prefactor for Yukawa term in fluid-fluid pair interactions");

    TwoDArray<double> YukawaKWF_Array( (Fluid_List->get<int>("F1_Ncomp"))*(Surface_List->get<int>("S3: Number of surface types")),1.0);
    PotentialsWF_List->set("WF5: ExpDecayParam_wf", YukawaKWF_Array, "alpha - exponential parameter in Yukawa term of fluid-fluid pair interactions");

    return;
}
/***********************************************************************************************************************************/

void dft_GUI_potentialsWF_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Surface_List, 
                         Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsWF_List )
{
  int i,j;
    /********************/
    /* setup parameters */
    /********************/
  /* can only set these parameters for manual mixing */
  if (Mix_type==1){
    TwoDArray<double> SigmaWF_Array(Ncomp,Nwall_type);
    for (i=0;i<Ncomp;i++) for (j=0;j<Nwall_type;j++) SigmaWF_Array[i][j]=Sigma_wf[i][j];
    PotentialsWF_List->set("WF1: Sigma_wf", SigmaWF_Array, "Sigma - Characteristic diameter fluid-fluid pair interactions");

    TwoDArray<double> EpsWF_Array(Ncomp,Nwall_type);
    for (i=0;i<Ncomp;i++) for (j=0;j<Nwall_type;j++) EpsWF_Array[i][j]=Eps_wf[i][j];
    PotentialsWF_List->set("WF2: Eps_wf", EpsWF_Array, "Eps - Energy prefactor for fluid-fluid pair interactions");

    TwoDArray<double> CutWF_Array(Ncomp,Nwall_type);
    for (i=0;i<Ncomp;i++) for (j=0;j<Nwall_type;j++) CutWF_Array[i][j]=Cut_wf[i][j];
    PotentialsWF_List->set("WF3: Cut_wf", CutWF_Array, "rcut - Cutoff distance for fluid-fluid pair interactions");

    TwoDArray<double> EpsYukawaWF_Array(Ncomp,Nwall_type);
    for (i=0;i<Ncomp;i++) for (j=0;j<Nwall_type;j++) EpsYukawaWF_Array[i][j]=EpsYukawa_wf[i][j];
    PotentialsWF_List->set("WF4: EpsYukawa_wf", EpsYukawaWF_Array, "AYuk - Energy prefactor for Yukawa term in fluid-fluid pair interactions");

    TwoDArray<double> YukawaKWF_Array(Ncomp,Nwall_type);
    for (i=0;i<Ncomp;i++) for (j=0;j<Nwall_type;j++) YukawaKWF_Array[i][j]=YukawaK_wf[i][j];
    PotentialsWF_List->set("WF5: ExpDecayParam_wf", YukawaKWF_Array, "alpha - exponential parameter in Yukawa term of fluid-fluid pair interactions");
  }
  return;
}
/***********************************************************************************************************************************/

void dft_GUI_potentialsWF_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Surface_List, 
                         Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsWF_List )
{
    /************************/
    /* define dependencies  */
    /************************/

    /* Visual Dependencies ..... note that it isn't really possible to do the visual dependencies for the various potential types at this time
     * because the vext_type IDs are in string arrays.  So, we would need to have a dependency based on an array entry.  This isn't possible at this time.*/

   Dependency::ParameterEntryList PotWF2DColNumber_Dependents;
   PotWF2DColNumber_Dependents.insert(PotentialsWF_List->getEntryRCP("WF1: Sigma_wf"));
   PotWF2DColNumber_Dependents.insert(PotentialsWF_List->getEntryRCP("WF2: Eps_wf"));
   PotWF2DColNumber_Dependents.insert(PotentialsWF_List->getEntryRCP("WF3: Cut_wf"));
   PotWF2DColNumber_Dependents.insert(PotentialsWF_List->getEntryRCP("WF4: EpsYukawa_wf"));
   PotWF2DColNumber_Dependents.insert(PotentialsWF_List->getEntryRCP("WF5: ExpDecayParam_wf"));

   RCP<TwoDColDependency<int,double> > PotWF2DColNumber_Dep = rcp(
        new TwoDColDependency<int,double>(Surface_List->getEntryRCP("S3: Number of surface types"),PotWF2DColNumber_Dependents));

   Dependency::ParameterEntryList PotWF2DRowNumber_Dependents;
   PotWF2DRowNumber_Dependents.insert(PotentialsWF_List->getEntryRCP("WF1: Sigma_wf"));
   PotWF2DRowNumber_Dependents.insert(PotentialsWF_List->getEntryRCP("WF2: Eps_wf"));
   PotWF2DRowNumber_Dependents.insert(PotentialsWF_List->getEntryRCP("WF3: Cut_wf"));
   PotWF2DRowNumber_Dependents.insert(PotentialsWF_List->getEntryRCP("WF4: EpsYukawa_wf"));
   PotWF2DRowNumber_Dependents.insert(PotentialsWF_List->getEntryRCP("WF5: ExpDecayParam_wf"));

   RCP<TwoDRowDependency<int,double> > PotWF2DRowNumber_Dep = rcp(
        new TwoDRowDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"),PotWF2DRowNumber_Dependents));


  /************************************************/
  /* add the dependencies to the dependency sheet */
  /************************************************/
    depSheet_Tramonto->addDependency(PotWF2DColNumber_Dep);
    depSheet_Tramonto->addDependency(PotWF2DRowNumber_Dep);

  return;
}
/***********************************************************************************************************************************/
