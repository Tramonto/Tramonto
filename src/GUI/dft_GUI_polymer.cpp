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

   Dependency::ParameterParentMap PolymerDependents;
   PolymerDependents.insert(std::pair<std::string, RCP<ParameterList> >("P1: Npoly_comp", Polymer_List));
   PolymerDependents.insert(std::pair<std::string, RCP<ParameterList> >("P2: Nblock[ipol_comp]", Polymer_List));
   PolymerDependents.insert(std::pair<std::string, RCP<ParameterList> >("P3: Block[ipol_comp][iblock]", Polymer_List));
   PolymerDependents.insert(std::pair<std::string, RCP<ParameterList> >("P4: BlockType[ipol_comp][iblock]", Polymer_List));
   PolymerDependents.insert(std::pair<std::string, RCP<ParameterList> >("P5: Any Grafted polymers?", Polymer_List));
   PolymerDependents.insert(std::pair<std::string, RCP<ParameterList> >("P8: Polymer achitecture entry", Polymer_List));

   RCP<StringVisualDependency> PolyAll_Dep = rcp(new StringVisualDependency("F4_POLYMER_Functional", Functional_List, PolymerDependents, 
           tuple<std::string>("Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT","Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)","Polymer_JDC_iSAFT(comp)")));

   Dependency::ParameterParentMap GraftDependents;
   GraftDependents.insert(std::pair<std::string, RCP<ParameterList> >("P6: Grafted_wall_ID[ipol_comp]", Polymer_List));
   GraftDependents.insert(std::pair<std::string, RCP<ParameterList> >("P7: Grafted_wall_Density[ipol_comp]", Polymer_List));

   RCP<BoolVisualDependency> Graft_Dep = rcp(new BoolVisualDependency("P5: Any Grafted polymers?", Polymer_List, GraftDependents, true));

   Dependency::ParameterParentMap CMSDependents;
   CMSDependents.insert(std::pair<std::string, RCP<ParameterList> >("CMS1: N_CrFiles (CMS)", PolymerCMS_List));
   CMSDependents.insert(std::pair<std::string, RCP<ParameterList> >("CMS2: Cr_File_1 (CMS)", PolymerCMS_List));
   CMSDependents.insert(std::pair<std::string, RCP<ParameterList> >("CMS5: Cr Radius (CMS)", PolymerCMS_List));

   RCP<StringVisualDependency> CMS_Dep = rcp(new StringVisualDependency("F4_POLYMER_Functional", Functional_List, CMSDependents, tuple<std::string>("Polymer_CMS")));

   Dependency::ParameterParentMap CMSFileDependents;
   CMSFileDependents.insert(std::pair<std::string, RCP<ParameterList> >("CMS3: Cr_File_2 (CMS)", PolymerCMS_List));
   CMSFileDependents.insert(std::pair<std::string, RCP<ParameterList> >("CMS4: CrFac (CMS)", PolymerCMS_List));

   /* doesn't work */
/*   RCP<NumberVisualDependency<int> > CMSFile_Dep = rcp(new NumberVisualDependency<int>("CMS1: N_CrFiles (CMS)", PolymerCMS_List, CMSFileDependents, 
                                                       (funcCMS(PolymerCMS_List->get<int>("CMS1: N_CrFiles (CMS)") )  ));*/


   RCP<StringCondition> PolyFile_Con1 = rcp(new StringCondition("P8: Polymer achitecture entry", Polymer_List, tuple<std::string>("Read From File"),true));
   RCP<StringCondition> PolyFile_Con2 = rcp(new StringCondition("F4_POLYMER_Functional", Functional_List, 
           tuple<std::string>("Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT","Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)","Polymer_JDC_iSAFT(comp)"),true));
   Condition::ConditionList PolyFile_conList = tuple<RCP<Condition> >(PolyFile_Con1, PolyFile_Con2);
   RCP<AndCondition> PolyFile_andCon = rcp(new AndCondition(PolyFile_conList));

   RCP<ConditionVisualDependency> PolyFile_Dep = rcp(new ConditionVisualDependency(PolyFile_andCon, "P9: Polymer architecture filename", Polymer_List, true));


   Dependency::ParameterParentMap PolymerArrayLength_Dependents;
   PolymerDependents.insert(std::pair<std::string, RCP<ParameterList> >("P2: Nblock[ipol_comp]", Polymer_List));
   PolymerDependents.insert(std::pair<std::string, RCP<ParameterList> >("P3: Block[ipol_comp][iblock]", Polymer_List));
   PolymerDependents.insert(std::pair<std::string, RCP<ParameterList> >("P4: BlockType[ipol_comp][iblock]", Polymer_List));
   PolymerDependents.insert(std::pair<std::string, RCP<ParameterList> >("P5: Grafted_wall_ID[ipol_comp]", Polymer_List));
   PolymerDependents.insert(std::pair<std::string, RCP<ParameterList> >("P6: Grafted_wall_Density[ipol_comp]", Polymer_List));

   RCP<NumberArrayLengthDependency> PolyAllLength_Dep = rcp(
           new NumberArrayLengthDependency( "P1: Npoly_comp", Polymer_List, PolymerArrayLength_Dependents));

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
