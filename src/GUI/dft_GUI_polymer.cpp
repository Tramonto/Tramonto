/*#include "Teuchos_ParameterList.hpp"*/
using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

int funcCMS(int nCrfiles);

void dft_GUI_Polymer_set_defaults( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                  Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerCMS_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerArch_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerGraft_List)
{

     /**************************************/
     /* Define validators for this section.*/
     /**************************************/

   RCP<EnhancedNumberValidator<int> > DimValidator = rcp(new EnhancedNumberValidator<int>(0,2,1));
   RCP<EnhancedNumberValidator<int> > Ncrfile_Validator = rcp(new EnhancedNumberValidator<int>(0,2,1));
   RCP<StringValidator> PolyArch_Validator = rcp<StringValidator>(new StringValidator(tuple<std::string>(
                      "Read From File","Set up in GUI","Linear Chains - Automatic set-up","Symmetric Linear Chains - Automate set-up")));
   RCP<FileNameValidator> polyFile_Validator = rcp(new FileNameValidator);
   RCP<FileNameValidator> CrFile_Validator = rcp(new FileNameValidator);

     /*********************/
     /* set up parameters */
     /*********************/
   Polymer_List->set("P1: Npoly_comp", 1, "Number of polymer components");

   Array<int> Nblock_Array( (Polymer_List->get<int>("P1: Npoly_comp")),1);
   Polymer_List->set("P2: Nblock_per_polymer", Nblock_Array, "Number of different chemical blocks in each polymer component.");

   Polymer_List->set("P3: Max_Nblocks", 1, "Maximum Number of blocks among polymer components");

   TwoDArray<int> Block_ID_Array(Polymer_List->get<int>("P1: Npoly_comp"),Polymer_List->get<int>("P3: Max_Nblocks"),0);
   Polymer_List->set("P4: Nseg_perBlock", Block_ID_Array, "Number of segments for each block in each polymer component.  Note that some array entries (iblock>Nblock[ipol_comp]) are irrelevant");

   TwoDArray<int> SegTypePerBlock_Array(Polymer_List->get<int>("P1: Npoly_comp"),Polymer_List->get<int>("P3: Max_Nblocks"),0);
   Polymer_List->set("P5: SegType_perBlock", SegTypePerBlock_Array, "Enter the segment component IDs for each block. Note that some array entries (iblock>Nblock[ipol_comp]) are irrelevant.");

   Polymer_List->set("P6: Polymer achitecture entry", "Read From File", "Identify how the polymer architecture will be set up",PolyArch_Validator);

   Polymer_List->set("P7: Polymer architecture filename","","Enter the file name where the polymer architecture can be found.", polyFile_Validator);



   PolymerGraft_List->set("PG1: Grafted Polymers?", false, "Set to true if there are any polymers grafted to a surface.");

/*   Array<bool> GraftedPolymer_Array(Polymer_List->getEntryRCP("P1: Npoly_comp"),false);
   PolymerGraft_List->set("PG2: Grafted_polymer_TF", GraftedPolymer_Array, "Identify the polymers that will be grafted to a wall.");*/

   Array<int> GraftedWall_Array( (Polymer_List->get<int>("P1: Npoly_comp")),-1);
   PolymerGraft_List->set("PG3: Grafted_wall_ID[ipol_comp]", GraftedWall_Array, "Identify the wall to which each polymer is grafted.\n  The flag value of -1 indicates no grafting.");

   Array<double> GraftedDensity_Array( (Polymer_List->get<int>("P1: Npoly_comp")),0.0);
   PolymerGraft_List->set("PG4: Grafted_wall_Density[ipol_comp]", GraftedDensity_Array, "The desired density of grafted chains on the surface.");



   PolymerCMS_List->set("CMS1: N_CrFiles",1,"Number of direct correlation function files to be read",Ncrfile_Validator); 

   PolymerCMS_List->set("CMS2: Cr_File_1", "", "Enter filename for a file containing a direct correlation function",CrFile_Validator);

   PolymerCMS_List->set("CMS3: Cr_File_2", "", "Enter filename for a file containing a second direct correlation function (DCF).\n If two files are provided, Tramonto will take an average of the two files using a parameter Crfac to compute a new direct correlation function.i\n  Specifically the hybrid DCF will be Crfac*Cr_File_1+(1.0-Crfac)*Cr_File_2.",CrFile_Validator);

   PolymerCMS_List->set("CMS4: CrFac",1.0,"Factor used for mixing of two direct correlation functions from different files.\n Specifically the hybrid DCF will be Crfac*Cr_File_1+(1.0-Crfac)*Cr_File_2."); 

   TwoDArray<double> CrRad_Array(Fluid_List->get<int>("F1_Ncomp"),Fluid_List->get<int>("F1_Ncomp"),1.0);
   PolymerCMS_List->set("CMS5: Cr HSRadius",CrRad_Array,"Radius of direct correlation function. Enter a value for each component pair"); 


   PolymerArch_List->set("PA1: NSeg_tot",1,"Enter the total number of coarse-grained segments in the model system.");
   PolymerArch_List->set("PA2: Nbond_max",2,"Enter the maximum number of bonds associated with any segment in the system.");

   Array<int> NBondAll_Array(PolymerArch_List->get<int>("PA1: NSeg_tot"),2);
   PolymerArch_List->set("PA3: NBondsperSeg",NBondAll_Array,"Enter the number of bonds for each segment in the system.\n For end segments and isolated single sites enter Nbond=2."); 

   TwoDArray<int> BondAll_Array(PolymerArch_List->get<int>("PA1: NSeg_tot"),PolymerArch_List->get<int>("PA2: Nbond_max"),-1);
   PolymerArch_List->set("PA4: BondAll",BondAll_Array,"For each segment (i) in the problem, enter the segment ids (jseg) to which iseg is bonded.\n End segments should be flagged with a jseg id = -1.\n Single site components should have two entries, both set to -1.\n Note that since the array is set based on maximum number of bonds for any segment in the system, some entries may be irrelevant."); 

   TwoDArray<int> PolySym_Array(PolymerArch_List->get<int>("PA1: NSeg_tot"),PolymerArch_List->get<int>("PA2: Nbond_max"),-1);
   PolymerArch_List->set("PA5: PolySym",PolySym_Array,"If you have a symmetric polymer you can indicated symmetric bonds here.  Enter -1 in all fields to ignore symmetric bonds.  Otherwise, enter the bond IDs to enforce symmetry. This will allow for some speed up of the calculations.\n Test your indexing by setting all to -1. Note that some entries may be irrelevant."); 

  return;
}
/***************************************************************************/
int funcCMS(int nCrfiles)
{
   if (nCrfiles==2) return 2;
   else return -1;
}
/***************************************************************************/
void dft_GUI_Polymer_set_OldFormat( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                  Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerCMS_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerArch_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerGraft_List)
{
  bool bool_tmp;
  int i,j,k,max,counter;

     /**************************************/
     /* Define validators for this section.*/
     /**************************************/

   RCP<EnhancedNumberValidator<int> > DimValidator = rcp(new EnhancedNumberValidator<int>(0,2,1));
   RCP<EnhancedNumberValidator<int> > Ncrfile_Validator = rcp(new EnhancedNumberValidator<int>(0,2,1));
   RCP<StringValidator> PolyArch_Validator = rcp<StringValidator>(new StringValidator(tuple<std::string>(
                      "Read From File","Set up in GUI","Linear Chains - Automatic set-up","Symmetric Linear Chains - Automate set-up")));
   RCP<FileNameValidator> polyFile_Validator = rcp(new FileNameValidator);
   RCP<FileNameValidator> CrFile_Validator = rcp(new FileNameValidator);

     /*********************/
     /* set up parameters */
     /*********************/
   if (Type_poly!=NONE){
      Polymer_List->set("P1: Npoly_comp", Npol_comp, "Number of polymer components");

      Array<int> Nblock_Array(Nblock,Nblock+Npol_comp);
      Polymer_List->set("P2: Nblock_per_polymer", Nblock_Array, "Number of different chemical blocks in each polymer component.");

      max=0;
      for (i=0;i<Npol_comp;i++) if (Nblock[i]>max) max=Nblock[i];
      Polymer_List->set("P3: Max_Nblocks", max, "Maximum Number of blocks among polymer components");

      TwoDArray<int> Block_ID_Array(Npol_comp,max);
      for (i=0; i<Npol_comp;i++) for (j=0; j<max;j++) Block_ID_Array(i,j)=Nseg_per_block[i][j];
      Polymer_List->set("P4: Nseg_perBlock", Block_ID_Array, "Number of segments for each block in each polymer component.  Note that some array entries (iblock>Nblock[ipol_comp]) are irrelevant");

      TwoDArray<int> SegTypePerBlock_Array(Npol_comp,max);
      for (i=0; i<Npol_comp;i++) for (j=0; j<max;j++) SegTypePerBlock_Array(i,j)=SegType_per_block[i][j];
      Polymer_List->set("P5: SegType_perBlock", SegTypePerBlock_Array, "Enter the segment component IDs for each block. Note that some array entries (iblock>Nblock[ipol_comp]) are irrelevant.");

      if (Type_poly_arch==POLY_ARCH_FILE)    Polymer_List->set("P6: Polymer achitecture entry", "Read From File", "Identify how the polymer architecture will be set up",PolyArch_Validator);
      else if (Type_poly_arch==LIN_POLY)     Polymer_List->set("P6: Polymer achitecture entry", "Linear Chains - Automatic set-up", "Identify how the polymer architecture will be set up",PolyArch_Validator);
      else if (Type_poly_arch==LIN_POLY_SYM) Polymer_List->set("P6: Polymer achitecture entry", "Symmetric Linear Chains - Automate set-up", "Identify how the polymer architecture will be set up",PolyArch_Validator);
      else                                   Polymer_List->set("P6: Polymer achitecture entry", "Set up in GUI", "Identify how the polymer architecture will be set up",PolyArch_Validator);

      if (Type_poly_arch==POLY_ARCH_FILE)
               Polymer_List->set("P7: Polymer architecture filename",(string)Poly_file_name,"Enter the file name where the polymer architecture can be found.", polyFile_Validator);
      else     Polymer_List->set("P7: Polymer architecture filename","","Enter the file name where the polymer architecture can be found.", polyFile_Validator);

      bool_tmp=false;
      for (i=0;i<Npol_comp;i++) if (Grafted[i]==TRUE) bool_tmp=true;
      PolymerGraft_List->set("PG1: Grafted Polymers?", bool_tmp, "Set to true if there are any polymers grafted to a surface.");

      Array<int> GraftedWall_Array(Graft_wall,Graft_wall+Npol_comp);
      PolymerGraft_List->set("PG3: Grafted_wall_ID[ipol_comp]", GraftedWall_Array, "Identify the wall to which each polymer is grafted.\n  The flag value of -1 indicates no grafting.");

      Array<double> GraftedDensity_Array(Rho_g,Rho_g+Npol_comp);
      PolymerGraft_List->set("PG4: Grafted_wall_Density[ipol_comp]", GraftedDensity_Array, "The desired density of grafted chains on the surface.");

      if (Type_poly == CMS || Type_poly==SCFT){
         PolymerCMS_List->set("CMS1: N_CrFiles",Ncr_files,"Number of direct correlation function files to be read",Ncrfile_Validator); 

         PolymerCMS_List->set("CMS2: Cr_File_1", (string)Cr_file, "Enter filename for a file containing a direct correlation function",CrFile_Validator);

         if (Ncr_files>1) PolymerCMS_List->set("CMS3: Cr_File_2", (string)Cr_file2, "Enter filename for a file containing a second direct correlation function (DCF).\n If two files are provided, Tramonto will take an average of the two files using a parameter Crfac to compute a new direct correlation function.i\n  Specifically the hybrid DCF will be Crfac*Cr_File_1+(1.0-Crfac)*Cr_File_2.",CrFile_Validator);
         else              PolymerCMS_List->set("CMS3: Cr_File_2", "", "Enter filename for a file containing a second direct correlation function (DCF).\n If two files are provided, Tramonto will take an average of the two files using a parameter Crfac to compute a new direct correlation function.i\n  Specifically the hybrid DCF will be Crfac*Cr_File_1+(1.0-Crfac)*Cr_File_2.",CrFile_Validator);

         PolymerCMS_List->set("CMS4: CrFac",Crfac,"Factor used for mixing of two direct correlation functions from different files.\n Specifically the hybrid DCF will be Crfac*Cr_File_1+(1.0-Crfac)*Cr_File_2."); 

         TwoDArray<double> CrRad_Array(Ncomp,Ncomp);
         for (i=0;i<Ncomp;i++) for (j=0;j<Ncomp;j++) CrRad_Array[i][j]=Cr_rad_hs[i][j];
         PolymerCMS_List->set("CMS5: Cr HSRadius",CrRad_Array,"Radius of direct correlation function. Enter a value for each component pair"); 

      }


      for (i=0; i<Npol_comp; i++) Nseg_tot += Nmer[i];

      PolymerArch_List->set("PA1: NSeg_tot",Nseg_tot,"Enter the total number of coarse-grained segments in the model system.");

      PolymerArch_List->set("PA2: Nbond_max",Nbond_max,"Enter the maximum number of bonds associated with any segment in the system.");

      Array<int> NBondAll_Array(Nseg_tot);
      counter=0;
      for (i=0; i<Npol_comp;i++)
         for (j=0; j<Nmer[i];j++){
            NBondAll_Array[counter]=Nbond[i][j];
            counter++;
      }
      PolymerArch_List->set("PA3: NBondsperSeg",NBondAll_Array,"Enter the number of bonds for each segment in the system.\n For end segments and isolated single sites enter Nbond=2."); 

      TwoDArray<int> BondAll_Array(Nseg_tot,Nbond_max);
      counter=0;
      for (i=0; i<Npol_comp;i++)
         for (j=0; j<Nmer[i];j++){
            for (k=0; k<Nbond[i][j];k++){
            BondAll_Array[counter][k]=Bonds[i][j][k];
            }
            counter++;
      }
      PolymerArch_List->set("PA4: BondAll",BondAll_Array,"For each segment (i) in the problem, enter the segment ids (jseg) to which iseg is bonded.\n End segments should be flagged with a jseg id = -1.\n Single site components should have two entries, both set to -1.\n Note that since the array is set based on maximum number of bonds for any segment in the system, some entries may be irrelevant."); 

      TwoDArray<int> PolySym_Array(Nseg_tot,Nbond_max);

      counter=0;
      for (i=0; i<Npol_comp;i++)
         for (j=0; j<Nmer[i];j++){
            for (k=0; k<Nbond[i][j];k++){
            PolySym_Array[counter][k]=pol_sym_tmp[i][j][k];
            }
            counter++;
      }
      PolymerArch_List->set("PA5: PolySym",PolySym_Array,"If you have a symmetric polymer you can indicated symmetric bonds here.  Enter -1 in all fields to ignore symmetric bonds.  Otherwise, enter the bond IDs to enforce symmetry. This will allow for some speed up of the calculations.\n Test your indexing by setting all to -1. Note that some entries may be irrelevant."); 

  }
  return;
}
/***************************************************************************/
void dft_GUI_Polymer_dependencies( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                  Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerCMS_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerArch_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerGraft_List)
{
        /************************/
        /* set up dependencies */
        /************************/

   Dependency::ParameterEntryList PolymerDependents;
   PolymerDependents.insert(Polymer_List->getEntryRCP("P1: Npoly_comp"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P2: Nblock_per_polymer"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P3: Max_Nblocks"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P4: Nseg_perBlock"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P5: SegType_perBlock"));
   PolymerDependents.insert(Polymer_List->getEntryRCP("P6: Polymer achitecture entry"));
   PolymerDependents.insert(PolymerGraft_List->getEntryRCP("PG1: Grafted Polymers?"));
   PolymerDependents.insert(PolymerGraft_List->getEntryRCP("PG3: Grafted_wall_ID[ipol_comp]"));
   PolymerDependents.insert(PolymerGraft_List->getEntryRCP("PG4: Grafted_wall_Density[ipol_comp]"));

   RCP<StringVisualDependency> PolyAll_Dep = rcp(new StringVisualDependency(Functional_List->getEntryRCP("F4_POLYMER_Functional"), PolymerDependents, 
           tuple<std::string>("Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT","Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)","Polymer_JDC_iSAFT(comp)")));

   Dependency::ParameterEntryList GraftDependents;
/*   GraftDependents.insert(PolymerGraft_List->getEntryRCP("PG2: Grafted_polymer_TF"));*/
   GraftDependents.insert(PolymerGraft_List->getEntryRCP("PG3: Grafted_wall_ID[ipol_comp]"));
   GraftDependents.insert(PolymerGraft_List->getEntryRCP("PG4: Grafted_wall_Density[ipol_comp]"));

   RCP<BoolVisualDependency> Graft_Dep = rcp(
        new BoolVisualDependency(PolymerGraft_List->getEntryRCP("PG1: Grafted Polymers?"), GraftDependents, true));

   Dependency::ParameterEntryList ArchVisDependents;
   ArchVisDependents.insert(PolymerArch_List->getEntryRCP("PA1: NSeg_tot"));
   ArchVisDependents.insert(PolymerArch_List->getEntryRCP("PA2: Nbond_max"));
   ArchVisDependents.insert(PolymerArch_List->getEntryRCP("PA3: NBondsperSeg"));
   ArchVisDependents.insert(PolymerArch_List->getEntryRCP("PA4: BondAll"));
   ArchVisDependents.insert(PolymerArch_List->getEntryRCP("PA5: PolySym"));
   RCP<StringVisualDependency> ArchVis_Dep = rcp(new StringVisualDependency(Polymer_List->getEntryRCP("P6: Polymer achitecture entry"), ArchVisDependents, tuple<std::string>("Set up in GUI")));

   Dependency::ParameterEntryList CMSDependents;
   CMSDependents.insert(PolymerCMS_List->getEntryRCP("CMS1: N_CrFiles"));
   CMSDependents.insert(PolymerCMS_List->getEntryRCP("CMS2: Cr_File_1"));
   CMSDependents.insert(PolymerCMS_List->getEntryRCP("CMS3: Cr_File_2"));
   CMSDependents.insert(PolymerCMS_List->getEntryRCP("CMS4: CrFac"));
   CMSDependents.insert(PolymerCMS_List->getEntryRCP("CMS5: Cr HSRadius"));
   RCP<StringVisualDependency> CMS_Dep = rcp(new StringVisualDependency(Functional_List->getEntryRCP("F4_POLYMER_Functional"), CMSDependents, tuple<std::string>("Polymer_CMS")));

   Dependency::ParameterEntryList CMSFileDependents;
   CMSFileDependents.insert(PolymerCMS_List->getEntryRCP("CMS3: Cr_File_2"));
   CMSFileDependents.insert(PolymerCMS_List->getEntryRCP("CMS4: CrFac"));
   /* doesn't work */
   RCP<NumberVisualDependency<int> > CMSFile_Dep = rcp(new NumberVisualDependency<int>(PolymerCMS_List->getEntryRCP("CMS1: N_CrFiles"), CMSFileDependents, 
                                                       (funcCMS(PolymerCMS_List->get<int>("CMS1: N_CrFiles") )  )));


   RCP<StringCondition> PolyFile_Con1 = rcp(new StringCondition(Polymer_List->getEntryRCP("P6: Polymer achitecture entry"), tuple<std::string>("Read From File")));
   RCP<StringCondition> PolyFile_Con2 = rcp(new StringCondition(Functional_List->getEntryRCP("F4_POLYMER_Functional"), 
           tuple<std::string>("Polymer_CMS","Polymer_CMS_SCFT","Polymer_TC_iSAFT","Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)","Polymer_JDC_iSAFT(comp)")));
   Condition::ConstConditionList PolyFile_conList = tuple<RCP<const Condition> >(PolyFile_Con1, PolyFile_Con2);
   RCP<AndCondition> PolyFile_andCon = rcp(new AndCondition(PolyFile_conList));

   RCP<ConditionVisualDependency> PolyFile_Dep = rcp(
       new ConditionVisualDependency(PolyFile_andCon,Polymer_List->getEntryRCP("P7: Polymer architecture filename"), true));

/*   Dependency::ParameterEntryList PolymerBoolArrayLength_Deps;
   PolymerBoolArrayLength_Deps.insert(Polymer_List->getEntryRCP("PG2: Grafted_polymer_TF"));
   RCP<NumberArrayLengthDependency<int,bool> > PolyBoolLength_Dep = rcp(
           new NumberArrayLengthDependency<int,bool>(Polymer_List->getEntryRCP("P1: Npoly_comp"), PolymerBoolArrayLength_Deps));*/

   Dependency::ParameterEntryList PolymerIntArrayLength_Deps;
   PolymerIntArrayLength_Deps.insert(Polymer_List->getEntryRCP("P2: Nblock_per_polymer"));
   PolymerIntArrayLength_Deps.insert(PolymerGraft_List->getEntryRCP("PG3: Grafted_wall_ID[ipol_comp]"));
   RCP<NumberArrayLengthDependency<int,int> > PolyIntLength_Dep = rcp(
           new NumberArrayLengthDependency<int,int>(Polymer_List->getEntryRCP("P1: Npoly_comp"), PolymerIntArrayLength_Deps));

   Dependency::ParameterEntryList PolymerDoubleArrayLength_Deps;
   PolymerDoubleArrayLength_Deps.insert(PolymerGraft_List->getEntryRCP("PG4: Grafted_wall_Density[ipol_comp]"));
   RCP<NumberArrayLengthDependency<int,double> > PolyDoubleLength_Dep = rcp(
           new NumberArrayLengthDependency<int,double>(Polymer_List->getEntryRCP("P1: Npoly_comp"), PolymerDoubleArrayLength_Deps));

   Dependency::ParameterEntryList Polymer2DRowNumber_Dependents;
      Polymer2DRowNumber_Dependents.insert(Polymer_List->getEntryRCP("P4: Nseg_perBlock"));
      Polymer2DRowNumber_Dependents.insert(Polymer_List->getEntryRCP("P5: SegType_perBlock"));
      RCP<TwoDRowDependency<int,int> > Polymer2DRowNumber_Dep = rcp(
                new TwoDRowDependency<int,int>(Polymer_List->getEntryRCP("P1: Npoly_comp"),Polymer2DRowNumber_Dependents));

   Dependency::ParameterEntryList Polymer2DColNumber_Dependents;
      Polymer2DColNumber_Dependents.insert(Polymer_List->getEntryRCP("P4: Nseg_perBlock"));
      Polymer2DColNumber_Dependents.insert(Polymer_List->getEntryRCP("P5: SegType_perBlock"));
      RCP<TwoDColDependency<int,int> > Polymer2DColNumber_Dep = rcp(
                new TwoDColDependency<int,int>(Polymer_List->getEntryRCP("P3: Max_Nblocks"),Polymer2DColNumber_Dependents));

   Dependency::ParameterEntryList PolymerCMS2DRowNumber_Deps;
      PolymerCMS2DRowNumber_Deps.insert(PolymerCMS_List->getEntryRCP("CMS5: Cr HSRadius"));
      RCP<TwoDRowDependency<int,double> > PolymerCMS2DRowNumber_Dep = rcp(
                new TwoDRowDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"),PolymerCMS2DRowNumber_Deps));

   Dependency::ParameterEntryList PolymerCMS2DColNumber_Deps;
      PolymerCMS2DColNumber_Deps.insert(PolymerCMS_List->getEntryRCP("CMS5: Cr HSRadius"));
      RCP<TwoDColDependency<int,double> > PolymerCMS2DColNumber_Dep = rcp(
                new TwoDColDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"),PolymerCMS2DColNumber_Deps));

   Dependency::ParameterEntryList PolymerArchArrayLength_Deps;
   PolymerArchArrayLength_Deps.insert(PolymerArch_List->getEntryRCP("PA3: NBondsperSeg"));
   RCP<NumberArrayLengthDependency<int,int> > PolyArchLength_Dep = rcp(
           new NumberArrayLengthDependency<int,int>(PolymerArch_List->getEntryRCP("PA1: NSeg_tot"), PolymerArchArrayLength_Deps));

   Dependency::ParameterEntryList PolymerArch2DRowNumber_Deps;
      PolymerArch2DRowNumber_Deps.insert(PolymerArch_List->getEntryRCP("PA4: BondAll"));
      PolymerArch2DRowNumber_Deps.insert(PolymerArch_List->getEntryRCP("PA5: PolySym"));
      RCP<TwoDRowDependency<int,int> > PolymerArch2DRowNumber_Dep = rcp(
                new TwoDRowDependency<int,int>(PolymerArch_List->getEntryRCP("PA1: NSeg_tot"),PolymerArch2DRowNumber_Deps));

   Dependency::ParameterEntryList PolymerArch2DColNumber_Deps;
      PolymerArch2DColNumber_Deps.insert(PolymerArch_List->getEntryRCP("PA4: BondAll"));
      PolymerArch2DColNumber_Deps.insert(PolymerArch_List->getEntryRCP("PA5: PolySym"));
      RCP<TwoDColDependency<int,int> > PolymerArch2DColNumber_Dep = rcp(
                new TwoDColDependency<int,int>(PolymerArch_List->getEntryRCP("PA2: Nbond_max"),PolymerArch2DColNumber_Deps));

      /*****************************************/
      /* add the dependencies for this section.*/
      /*****************************************/

   depSheet_Tramonto->addDependency(PolyAll_Dep);
   depSheet_Tramonto->addDependency(Graft_Dep);
   depSheet_Tramonto->addDependency(CMS_Dep);
   depSheet_Tramonto->addDependency(CMSFile_Dep);
   depSheet_Tramonto->addDependency(ArchVis_Dep);
   depSheet_Tramonto->addDependency(PolyFile_Dep);
   /*depSheet_Tramonto->addDependency(PolyBoolLength_Dep);*/
   depSheet_Tramonto->addDependency(PolyIntLength_Dep);
   depSheet_Tramonto->addDependency(PolyDoubleLength_Dep);
   depSheet_Tramonto->addDependency(Polymer2DRowNumber_Dep);
   depSheet_Tramonto->addDependency(Polymer2DColNumber_Dep);
   depSheet_Tramonto->addDependency(PolymerCMS2DRowNumber_Dep);
   depSheet_Tramonto->addDependency(PolymerCMS2DColNumber_Dep);
   depSheet_Tramonto->addDependency(PolyArchLength_Dep);
   depSheet_Tramonto->addDependency(PolymerArch2DRowNumber_Dep);
   depSheet_Tramonto->addDependency(PolymerArch2DColNumber_Dep);

   return;
}
