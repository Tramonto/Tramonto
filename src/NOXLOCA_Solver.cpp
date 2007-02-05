//@HEADER
// ********************************************************************
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER

#include "LOCA.H"
#include "NOXLOCA_Tramonto_Group.hpp"
#include "NOXLOCA_Tramonto_PTGroup.hpp"

extern "C" {
#include "loca_const.h"

void NOXLOCA_Solver(double **xBox, double **xOwned, double **x2Owned)
{
  try {
  
    // Create a solution vector of NOX's known type
    NOXLOCA::Tramonto::Vector xTV(Nunk_per_node, Nnodes_per_proc, xOwned);

    // Create parameter list
    Teuchos::RefCountPtr<Teuchos::ParameterList> paramList = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    
    if (Loca.method == 2)  stepperList.set("Continuation Method", "Arc Length");
    else  stepperList.set("Continuation Method", "Natural");


    // Set up the problem interface
    LOCA::ParameterVector p;
    double init_bif_param = -1.0;

    // Newton solve -- create dummy parameter that is ignored
    if (Loca.method < 0) {
      p.addParameter("NoParam", 0.0);
      stepperList.set("Continuation Parameter", "NoParam");
      stepperList.set("Initial Value", 0.0);
      stepperList.set("Max Steps", 0);
    }
    // Continuation run
    else {
      double init_param =  get_init_param_value(Loca.cont_type1);
      p.addParameter("ConParam", init_param);
      stepperList.set("Continuation Parameter", "ConParam");
      stepperList.set("Initial Value", init_param);
      stepperList.set("Max Steps", Loca.num_steps);

      // Truning point (spinodal) run or phase transition tracking run
      if (Loca.method == 3 || Loca.method == 4) {
        init_bif_param = get_init_param_value(Loca.cont_type2);
        p.addParameter("BifParam", init_bif_param);
      }

if (Loca.method == 4) { cout << "NO PHASE TRANSITION ALG just double-solving so far!" << endl; }
    }

    stepperList.set("Max Value", 1.0e8);
    stepperList.set("Min Value", -1.0e8);
    stepperList.set("Max Nonlinear Iterations", Max_Newton_iter);
    stepperList.set("Enable Arc Length Scaling", true);
    stepperList.set("Goal Arc Length Parameter Contribution", 0.5);
    stepperList.set("Max Arc Length Parameter Contribution", 0.7);
    stepperList.set("Initial Scale Factor", 1.0);
    stepperList.set("Min Scale Factor", 1.0e-8);
    stepperList.set("Enable Tangent Factor Step Size Scaling",true);
    stepperList.set("Min Tangent Factor", -1.0);
    stepperList.set("Tangent Factor Exponent",1.0);
    stepperList.set("Compute Eigenvalues",false);
    stepperList.set("Bordered Solver Method", "Bordering");

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    if (Loca.method == 3) {
       bifurcationList.set("Type", "Turning Point");
       bifurcationList.set("Bifurcation Parameter", "BifParam");
       Teuchos::RefCountPtr<NOX::Abstract::Vector> nullVec =
         Teuchos::rcp(new NOXLOCA::Tramonto::Vector(xTV));
       nullVec->init(1.0);
       bifurcationList.set("Length Normalization Vector", nullVec);
       bifurcationList.set("Initial Null Vector", nullVec);
       bifurcationList.set("Update Null Vectors Every Continuation Step", true);
   }
    else bifurcationList.set("Type", "None");

    // Create predictor sublist
    Teuchos::ParameterList& predictorList = locaParamsList.sublist("Predictor");
    if (Loca.method == 0 )
      predictorList.set("Method", "Constant");
    else if (Loca.method == 4  || Loca.method == 3)
      predictorList.set("Method", "Secant");
    else
      predictorList.set("Method", "Tangent");
    //predictorList.set("Method", "Secant");

    //Teuchos::ParameterList& firstStepPredictorList = 
    //  predictorList.sublist("First Step Predictor");
    //firstStepPredictorList.set("Method", "Tangent");

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    //stepSizeList.set("Method", "Constant");
    stepSizeList.set("Method", "Adaptive");
    stepSizeList.set("Initial Step Size", Loca.step_size);
    //stepSizeList.set("Min Step Size", 1.0e-9);
    //stepSizeList.set("Max Step Size", 1.0e5);
    stepSizeList.set("Aggressiveness", Loca.aggr); // for adaptive
    stepSizeList.set("Failed Step Reduction Factor", 0.5);
    stepSizeList.set("Successful Step Increase Factor", 1.26); // for constant

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    nlParams.set("Nonlinear Solver", "Line Search Based");
      Teuchos::ParameterList& lineSearchParameters = nlParams.sublist("Line Search");
        lineSearchParameters.set("Method","Backtrack"); //Cuts step for NANs
        Teuchos::ParameterList& backtrackParameters = 
            lineSearchParameters.sublist("Backtrack");
        backtrackParameters.set("Minimum Step", 0.005);

    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("Output Information", 
			       NOX::Utils::Details +
			       NOX::Utils::OuterIteration + 
			       NOX::Utils::InnerIteration + 
			       NOX::Utils::Warning + 
			       NOX::Utils::StepperIteration +
			       NOX::Utils::StepperDetails);


    // Create global data object
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList);

    // Create the Group -- this is the main guy
    Teuchos::RefCountPtr<NOXLOCA::Tramonto::Group> grp = 
      Teuchos::rcp(new NOXLOCA::Tramonto::Group(globalData, xTV, xBox, p, paramList));

    Teuchos::RefCountPtr<LOCA::Abstract::Group> solveGrp; // Pointer, set in if block

    if (Loca.method == 4) {
      NOXLOCA::Tramonto::Vector x2TV(Nunk_per_node, Nnodes_per_proc, x2Owned);

      NOXLOCA::Tramonto::PTVector ptVec(xTV, x2TV, init_bif_param);
      Teuchos::RefCountPtr<NOXLOCA::Tramonto::PTGroup> ptgrp = 
        Teuchos::rcp(new NOXLOCA::Tramonto::PTGroup(globalData, grp, ptVec));
 
      solveGrp = ptgrp;
    }
    else {
      solveGrp = grp; // All execpt phase transitions
    }
    
    // Set up the status tests
    double tol=1.0e-9; if (Loca.method==3) tol=1.0e-6;

    Teuchos::RefCountPtr<NOX::StatusTest::NormF> normF = 
      Teuchos::rcp(new NOX::StatusTest::NormF(tol));
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxIters = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(Max_Newton_iter));
    Teuchos::RefCountPtr<NOX::StatusTest::Generic> comboOR = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					      normF, 
					      maxIters));

    Teuchos::RefCountPtr<Teuchos::ParameterList> nlParamsRCP =
	    Teuchos::rcp(&nlParams, false);

    // Create the stepper  
    LOCA::Stepper stepper(globalData, solveGrp, comboOR, paramList);

    // Perform continuation run
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    // Check for convergence
    if (status != LOCA::Abstract::Iterator::Finished) {
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
	globalData->locaUtils->out() 
	  << "Stepper failed to converge!" << std::endl;
    }

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::Parameters)) {
      globalData->locaUtils->out() 
	<< std::endl << "Final Parameters" << std::endl
	<< "****************" << std::endl;
      paramList->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    destroyGlobalData(globalData);
  }

  catch (std::exception& e) {
    cout << e.what() << endl;
  }
  catch (const char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

}

}
