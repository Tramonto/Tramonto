//using namespace std;
#include <iostream>
#include "dft_globals_const.h" 
#include "dft_GUI.h"
#include "dft_GUI.hpp"

using namespace Teuchos;
using namespace Optika;

extern "C" void dft_OptikaGUI_control()
{
  int i,n,m;
  bool read_xml_file=false, read_input_old_format=false,start_GUI=true;
  string runpath, inputOLD_File,inputXML_File, inputGUI_File, output_File, outpath;


  /* set up a parameter list for a very simple GUI to select the type of run we are interested in */
  RCP<ParameterList> RunType_List = RCP<ParameterList>((new ParameterList("Root RunType List")));
  RCP<DependencySheet> depSheet_RunType = rcp(new DependencySheet());

  RCP<StringValidator> RunTypeValidator = rcp(
        new StringValidator(tuple<std::string>(
           "GUI - Defaults", 
           "GUI - Old Format File Parameters",
           "GUI - XML File Parameters",
           "No GUI / Old Format File",
           "No GUI / XML File")));
 
  const bool useFileDialogForDir = false;

  RCP<FileNameValidator> inputXMLFileVal = rcp (new FileNameValidator);
  RCP<FileNameValidator> inputOLDFileVal = rcp (new FileNameValidator);
  RCP<FileNameValidator> inputGUIFileVal = rcp (new FileNameValidator);

  inputXMLFileVal->setFileMustExist(true);
  inputXMLFileVal->setFileEmptyNameOK(true);
  inputOLDFileVal->setFileMustExist(true);
  inputOLDFileVal->setFileEmptyNameOK(true);

  if(useFileDialogForDir)
    inputGUIFileVal->setFileMustExist(true);

  inputGUIFileVal->setFileEmptyNameOK(true);

  RCP<FileNameValidator> outputFileVal = rcp (new FileNameValidator);

  RunType_List->set("R1: Run Type","GUI - Defaults","Indicate how you would like to run Tramonto for this job.",RunTypeValidator);

  RunType_List->set("R2: Input File (OLD format)","","select the OLD format input file. This choice also sets run directory.",inputOLDFileVal);
  RunType_List->set("R2: Input File (XML)","","select the XML input file. This choice also sets run directory.",inputXMLFileVal);

  if(useFileDialogForDir)
    RunType_List->set("R2: Input Directory (Set any File)","","select any file in desired directory.  The runpath will be extracted from the filename.",inputGUIFileVal);
  else
    RunType_List->set("R2: Input Directory","","select any input directory.",inputGUIFileVal);

  if(useFileDialogForDir)
    RunType_List->set("R3: Output Directory (Set any File)","","set any file in a desired output directory.  The output path will be extracted from the filename.",outputFileVal);
  else
    RunType_List->set("R3: Output Directory","","set any output directory.",outputFileVal);

   RCP<StringVisualDependency> RunType_Dep1 = rcp(new StringVisualDependency(
          RunType_List->getEntryRCP("R1: Run Type"),
          useFileDialogForDir ?
            RunType_List->getEntryRCP("R2: Input Directory (Set any File)") :
            RunType_List->getEntryRCP("R2: Input Directory"),
          tuple<std::string>("GUI - Defaults")));

   RCP<StringVisualDependency> RunType_Dep2 = rcp(new StringVisualDependency(
          RunType_List->getEntryRCP("R1: Run Type"), RunType_List->getEntryRCP("R2: Input File (XML)"), 
          tuple<std::string>("GUI - XML File Parameters","No GUI / XML File")));

   RCP<StringVisualDependency> RunType_Dep3 = rcp(new StringVisualDependency(
          RunType_List->getEntryRCP("R1: Run Type"), RunType_List->getEntryRCP("R2: Input File (OLD format)"), 
          tuple<std::string>("GUI - Old Format File Parameters","No GUI / Old Format File")));

   depSheet_RunType->addDependency(RunType_Dep1);
   depSheet_RunType->addDependency(RunType_Dep2);
   depSheet_RunType->addDependency(RunType_Dep3);
  
 /* Greate the GUI window for run startup for Tramonto  */
  Optika::getInput(RunType_List,depSheet_RunType);


  if (RunType_List->get<string>("R1: Run Type")=="GUI - Defaults"){
    if(useFileDialogForDir) {
      output_File=RunType_List->get<string>("R3: Output Directory (Set any File)");
      string::size_type m=output_File.find_last_of("/");
      outpath=output_File.substr(0,m);

      inputGUI_File=RunType_List->get<string>("R2: Select Any File in Desired Directory");
      string::size_type n=inputGUI_File.find_last_of("/");
      runpath=inputGUI_File.substr(0,n);
    }
    else {
      outpath=RunType_List->get<string>("R3: Output Directory");
      runpath=RunType_List->get<string>("R2: Input Directory");
    }
  }
  else if (RunType_List->get<string>("R1: Run Type")=="GUI - Old Format File Parameters"){
    if(useFileDialogForDir) {
      output_File=RunType_List->get<string>("R3: Output Directory (Set any File)");
      string::size_type m=output_File.find_last_of("/");
      outpath=output_File.substr(0,m);
    }
    else {
      outpath=RunType_List->get<string>("R3: Output Directory");
    }

     inputOLD_File=RunType_List->get<string>("R2: Input File (OLD format)");
     string::size_type n=inputOLD_File.find_last_of("/");
     runpath=inputOLD_File.substr(0,n);
     read_input_old_format=true;
  }
  else if (RunType_List->get<string>("R1: Run Type")=="GUI - XML File Parameters"){
    if(useFileDialogForDir) {
      output_File=RunType_List->get<string>("R3: Output Directory (Set any File)");
      string::size_type m=output_File.find_last_of("/");
      outpath=output_File.substr(0,m);
    }
    else {
      outpath=RunType_List->get<string>("R3: Output Directory");
    }

     inputXML_File=RunType_List->get<string>("R2: Input File (XML)");
     string::size_type n=inputXML_File.find_last_of("/");
     runpath=inputXML_File.substr(0,n);
     read_xml_file=true;
  }
  else if (RunType_List->get<string>("R1: Run Type")=="No GUI / Old Format File"){
    if(useFileDialogForDir) {
      output_File=RunType_List->get<string>("R3: Output Directory (Set any File)");
      string::size_type m=output_File.find_last_of("/");
      outpath=output_File.substr(0,m);
    }
    else {
      outpath=RunType_List->get<string>("R3: Output Directory");
    }

     inputOLD_File=RunType_List->get<string>("R2: Input File (OLD format)");
     string::size_type n=inputOLD_File.find_last_of("/");
     runpath=inputOLD_File.substr(0,n);
     read_input_old_format=true;
     start_GUI=false;
  }
  else if (RunType_List->get<string>("R1: Run Type")=="No GUI / XML File"){
    if(useFileDialogForDir) {
      output_File=RunType_List->get<string>("R3: Output Directory (Set any File)");
      string::size_type m=output_File.find_last_of("/");
      outpath=output_File.substr(0,m);
    }
    else {
      outpath=RunType_List->get<string>("R3: Output Directory");
    }

     inputXML_File=RunType_List->get<string>("R2: Input File (XML)");
     string::size_type n=inputXML_File.find_last_of("/");
     runpath=inputXML_File.substr(0,n);
     read_xml_file=true;
     start_GUI=false;
  }
  runpath=runpath+"/";
  outpath=outpath+"/";


  if (read_input_old_format) Read_OLDInput_File=TRUE;
  else                       Read_OLDInput_File=FALSE;

  if (read_xml_file) Read_XMLInput_File=TRUE;
  else               Read_XMLInput_File=FALSE;

  if (start_GUI) Open_GUI=TRUE;  
  else           Open_GUI=FALSE;  

  Runpath=(char *)runpath.c_str();
  Outpath=(char *)outpath.c_str();
  InputXML_File=(char *)inputXML_File.c_str();
  InputOLD_File=(char *)inputOLD_File.c_str();

  return;
}

