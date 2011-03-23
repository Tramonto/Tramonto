/*#include "Teuchos_ParameterList.hpp"*/
using namespace std;
#include <iostream>
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

int funcCMS(int nCrfiles);

void dft_GUI_Polymer( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                  Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerCMS_List)
{
  int idim;


        /**************************************/
        /* Define validators for this section.*/
        /**************************************/

    RCP<EnhancedNumberValidator<int> > DimValidator = rcp(new EnhancedNumberValidator<int>(0,2,1));
    RCP<EnhancedNumberValidator<int> > Ncrfile_Validator = rcp(new EnhancedNumberValidator<int>(1,2,1));
    RCP<StringValidator> PolyArch_Validator = rcp<StringValidator>(new StringValidator(tuple<std::string>(
                       "Read From File","Set up in GUI","Linear Chains - Automatic set-up","Symmetric Linear Chains - Automate set-up")));
    RCP<FileNameValidator> polyFile_Validator = rcp(new FileNameValidator);
    RCP<FileNameValidator> CrFile_Validator = rcp(new FileNameValidator);



        /*********************/
        /* set up parameters */
        /*********************/

   Polymer_List->set("P1: Npoly_comp", 1, "Number of polymer components");

   Array<int> Nblock_Array( (Polymer_List->get<int>("P1: Npoly_comp")),1);
   Polymer_List->set("P2: Nblock[ipol_comp]", Nblock_Array, "Number of different chemical blocks in each polymer component.");

   Array<int> Block_Array( (Polymer_List->get<int>("P1: Npoly_comp")),0);
   Polymer_List->set("P3: Block[ipol_comp][iblock]", Block_Array, "Block ids in each polymer component.");

   Array<int> BlockType_Array( (Polymer_List->get<int>("P1: Npoly_comp")),0);
   Polymer_List->set("P4: BlockType[ipol_comp][iblock]", BlockType_Array, "Block types should correspond to segment component ids.");

   Polymer_List->set("P5: Any Grafted polymers?", false, "Set to true if there are any polymers grafted to a surface.");

   Array<int> GraftedWall_Array( (Polymer_List->get<int>("P1: Npoly_comp")),-1);
   Polymer_List->set("P6: Grafted_wall_ID[ipol_comp]", GraftedWall_Array, "Identify the wall to which each polymer is grafted.\n  The flag value of -1 indicates no grafting.");

   Array<double> GraftedDensity_Array( (Polymer_List->get<int>("P1: Npoly_comp")),0.0);
   Polymer_List->set("P7: Grafted_wall_Density[ipol_comp]", GraftedDensity_Array, "The desired density of grafted chains on the surface.");

   Polymer_List->set("P8: Polymer achitecture entry", "Read From File", "Identify how the polymer architecture will be set up",PolyArch_Validator);

   Polymer_List->set("P9: Polymer architecture filename","","Enter the file name where the polymer architecture can be found.", polyFile_Validator);
 
   PolymerCMS_List->set("CMS1: N_CrFiles (CMS)",1,"Number of direct correlation function files to be read",Ncrfile_Validator); 

   PolymerCMS_List->set("CMS2: Cr_File_1 (CMS)", "", "Enter filename for a file containing a direct correlation function",CrFile_Validator);

   PolymerCMS_List->set("CMS3: Cr_File_2 (CMS)", "", "Enter filename for a file containing a second direct correlation function (DCF).\n If two files are provided, Tramonto will take an average of the two files using a parameter Crfac to compute a new direct correlation function.i\n  Specifically the hybrid DCF will be Crfac*Cr_File_1+(1.0-Crfac)*Cr_File_2.",CrFile_Validator);

   PolymerCMS_List->set("CMS4: CrFac (CMS)",1.0,"Factor used for mixing of two direct correlation functions from different files.\n Specifically the hybrid DCF will be Crfac*Cr_File_1+(1.0-Crfac)*Cr_File_2."); 

   PolymerCMS_List->set("CMS5: Cr Radius (CMS)",1.0,"Radius of direct correlation function."); 

/* would be nice to have code here to set up chain architecture */
/* Need Nbond[ipol][iseg], Bond[ipol][iseg][ibond], PolSym[ipol][iseg][ibond] */
/* can make this Nbond[iseg_all], Bond[iseg_all][ibond], PolSym[iseg_all][ibond] */  
/* logically tricky to reduce to less than 2D array */


        /************************/
        /* set up dependencies */
        /************************/

   Dependency::ParameterEntryList PolymerDependents;
   PolymerDependents.insert(Polymer_List->getEntryRCP("P1: Npoly_comp"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P2: Nblock[ipol_comp]"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P3: Block[ipol_comp][iblock]"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P4: BlockType[ipol_comp][iblock]"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P5: Any Grafted polymers?"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P8: Polymer achitecture entry"));

   RCP<StringVisualDependency> PolyAll_Dep = rcp(
       new StringVisualDependency(
           Functional_List->getEntryRCP("F4_POLYMER_Functional"), 
           PolymerDependents, 
           tuple<std::string>("Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT",
             "Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)","Polymer_JDC_iSAFT(comp)"))
   );

   Dependency::ParameterEntryList GraftDependents;
   GraftDependents.insert(Polymer_List->getEntryRCP("P6: Grafted_wall_ID[ipol_comp]"));
   GraftDependents.insert(Polymer_List->getEntryRCP("P7: Grafted_wall_Density[ipol_comp]"));

   RCP<BoolVisualDependency> Graft_Dep = rcp(
       new BoolVisualDependency(
           Polymer_List->getEntryRCP("P5: Any Grafted polymers?"),
           GraftDependents, 
           true)
   );

   Dependency::ParameterEntryList CMSDependents;
   CMSDependents.insert(PolymerCMS_List->getEntryRCP("CMS1: N_CrFiles (CMS)"));
   CMSDependents.insert(PolymerCMS_List->getEntryRCP("CMS2: Cr_File_1 (CMS)"));
   CMSDependents.insert(PolymerCMS_List->getEntryRCP("CMS5: Cr Radius (CMS)"));

   RCP<StringVisualDependency> CMS_Dep = rcp(
       new StringVisualDependency(
           Functional_List->getEntryRCP("F4_POLYMER_Functional"), 
           CMSDependents, 
           tuple<std::string>("Polymer_CMS"))
   );

   Dependency::ParameterEntryList CMSFileDependents;
   CMSFileDependents.insert(PolymerCMS_List->getEntryRCP("CMS3: Cr_File_2 (CMS)"));
   CMSFileDependents.insert(PolymerCMS_List->getEntryRCP("CMS4: CrFac (CMS)"));

   /* doesn't work */
/*   RCP<NumberVisualDependency<int> > CMSFile_Dep = rcp(new NumberVisualDependency<int>("CMS1: N_CrFiles (CMS)", PolymerCMS_List, CMSFileDependents, 
                                                       (funcCMS(PolymerCMS_List->get<int>("CMS1: N_CrFiles (CMS)") )  ));*/


   RCP<StringCondition> PolyFile_Con1 = rcp(
       new StringCondition(
           Polymer_List->getEntryRCP("P8: Polymer achitecture entry"), 
           tuple<std::string>("Read From File"),
           true)
   );

   RCP<StringCondition> PolyFile_Con2 = rcp(
       new StringCondition(
           Functional_List->getEntryRCP("F4_POLYMER_Functional"), 
           tuple<std::string>("Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT",
               "Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)",
               "Polymer_JDC_iSAFT(comp)"),
           true)
   );

   Condition::ConstConditionList PolyFile_conList = tuple<RCP<const Condition> >(PolyFile_Con1, PolyFile_Con2);
   RCP<AndCondition> PolyFile_andCon = rcp(new AndCondition(PolyFile_conList));

   RCP<ConditionVisualDependency> PolyFile_Dep = rcp(
       new ConditionVisualDependency(
           PolyFile_andCon, 
           Polymer_List->getEntryRCP("P9: Polymer architecture filename"), 
           true)
   );


   Dependency::ParameterEntryList PolymerArrayLength_Dependents;
   PolymerDependents.insert(Polymer_List->getEntryRCP("P2: Nblock[ipol_comp]"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P3: Block[ipol_comp][iblock]"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P4: BlockType[ipol_comp][iblock]"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P5: Grafted_wall_ID[ipol_comp]"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P6: Grafted_wall_Density[ipol_comp]"));

   RCP<NumberArrayLengthDependency<int, int> > PolyAllLength_Dep = rcp(
       new NumberArrayLengthDependency<int, int> ( 
           Polymer_List->getEntryRCP("P1: Npoly_comp"), 
           PolymerArrayLength_Dependents)
   );

      /*****************************************/
      /* add the dependencies for this section.*/
      /*****************************************/

   depSheet_Tramonto->addDependency(PolyAll_Dep);
   depSheet_Tramonto->addDependency(Graft_Dep);
   depSheet_Tramonto->addDependency(CMS_Dep);
   depSheet_Tramonto->addDependency(PolyFile_Dep);
   depSheet_Tramonto->addDependency(PolyAllLength_Dep);

  return;
}
/***************************************************************************/
int funcCMS(int nCrfiles)
{
   if (nCrfiles==2) return 1;
   else return -1;
}
