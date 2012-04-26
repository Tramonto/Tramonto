using namespace std;
#include <iostream>
#include "dft_globals_const.h"
#include "dft_GUI.h"
#include "dft_GUI.hpp"
using namespace Teuchos;
using namespace Optika;

void dft_GUI_NumericalMethods_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Solver_List,
                         Teuchos::RCP<Teuchos::ParameterList> Coarsening_List,
                         Teuchos::RCP<Teuchos::ParameterList> LoadBalance_List,
                         Teuchos::RCP<Teuchos::ParameterList> PhysicsMethod_List,
                         Teuchos::RCP<Teuchos::ParameterList> NonlinearSolver_List,
                         Teuchos::RCP<Teuchos::ParameterList> LinearSolver_List)
{
  string str_coarse,str_zones,str_updatefrac,str_scalefac;

    /**************************************/
    /* STRINGS */
    /**************************************/
     str_coarse="Select a method for coarsening of the problem.\n options: bulk_zone: Set to true to replace DFT Euler-Lagrange equation with rho(r)=rho_b some distance from the surfaces.\n \t Poisson-Bolzmann zone: Set to true to replace DFT Euler-Lagrange equation with rho(r)=rho_b some distance from the surfaces.\n\t Residual coarsening: Coarsening of residual equations and Jacobian integrals in zones away from the surfaces.\n\t Jacobian coarsening: Coarsening of jacobian integrals only in zones away from surfaces.";

    str_zones= "Set the minimum distance from surface in each zone. \n The distances should be arranged from nearest to furthest where\n the nearest zone to the surface is the most refined zone,\n and the furthest zone from the surfaces is the most coarse zone.  The first entry should be 0.0";

    str_updatefrac= "Set a minimum fraction for solution updates if using the built in Newton solver.\n This allows code to perform fractional updates if convergence difficulties are identifed. Set to 1.0 for full Newton steps. \n For Picard solver, this fraction is a mixing parameter for old and new solutions.  If using the mixed Picard-Newton Solver enter the minimum Newton update fraction.\n The Picard mixing parameter will be set to the newton update parameter divided by 100.";

    str_scalefac="Enter the desired scaling factor for segment component on each polymer chain.  \nNote that there will be zero entries wherever there are no segments of a given type on a particular chain";
    /**************************************/
    /* VALIDATORS */
    /**************************************/
 
 
    RCP<StringValidator> NLSolverTypeValidator = rcp(
           new StringValidator(tuple<std::string>("Newton Built-In","Newton NOX", 
						  "Picard Built-In","Picard NOX",
						  "Picard/Newton Built-In","Picard/Newton NOX")));

    RCP<StringValidator> PhysicsScalingValidator = rcp(
           new StringValidator(tuple<std::string>("No Physics Scaling","Automatic Calculation","Manual Input")));

    RCP<StringValidator> CoarsenTypeValidator = rcp(
           new StringValidator(tuple<std::string>("none","Bulk zone","Poisson-Boltzmann zone",
                                  "residual and jacobian coarsening","jacobian only coarsening")));

    RCP<EnhancedNumberValidator<int> > Niter_Validator = rcp(new EnhancedNumberValidator<int>());

    RCP<EnhancedNumberValidator<double> > NLtolValidator = rcp(new EnhancedNumberValidator<double>(1.e-10,1.e-2,1.e-6,8));

    RCP<EnhancedNumberValidator<double> > NLupdate_Validator = rcp(new EnhancedNumberValidator<double>(0.001,1.0,0.01));

    RCP<StringValidator> LoadBalValidator = rcp(
           new StringValidator(tuple<std::string>("Linear Matrix Balance","Geometric Recursive Bisection","Weighted Recursive Bisection")));

    RCP<StringValidator> LoadBalWeightValidator = rcp(
           new StringValidator(tuple<std::string>("Set automatically","Make all weights equal")));

    RCP<StringValidator> LSolverTypeValidator = rcp(
           new StringValidator(tuple<std::string>("Schur Solver","General/AztecOO solver")));

    RCP<StringValidator> LSolverApproachValidator = rcp(
           new StringValidator(tuple<std::string>("GMRES","cg","tfqmr","cgs","bicgstab")));

    RCP<StringValidator> LScalingValidator = rcp(
           new StringValidator(tuple<std::string>("none","row_sum","jacobi","symrow_sum")));

    RCP<StringValidator> LPreconditionerValidator = rcp(
           new StringValidator(tuple<std::string>("none","ilu","ilut","Jacobi","symmetric Gauss-Seidel","LSpoly3")));

    RCP<EnhancedNumberValidator<int> > DimValidator = rcp(new EnhancedNumberValidator<int>(0,2,1));

    RCP<StringValidator> JacCoarsenTypeValidator = rcp(
         new StringValidator(tuple<std::string>("none",
        "Jacobian Coarsening identical to Residual Coarsening",
	"Jacobian: factor of 2 more coarse than resid in most refined zone",
	"Jacobian: factor of 2 more coarse than resid in all but most coarse zone",
	"Jacobian: use most coarse zone to define entire matrix",
	"Jacobian: use 2nd most coarse zone in all but most coarse region",
	"Jacobian: set Esize_jac for all matrix calculations")));


    /***************************************************************/
    /* set up parameters */
    /***************************************************************/
    LoadBalance_List->set("LB1: Load Balancing Approach", "Weighted Recursive Bisection", "Select a method for load balancing a parallel processing job", LoadBalValidator);
    LoadBalance_List->set("LB2: Load Balancing Weights", "Set automatically", "How to assign weights for weighted rcb method", LoadBalWeightValidator);

    Coarsening_List->set("C1: Coarsening Type", "none", str_coarse, CoarsenTypeValidator);

    Coarsening_List->set("C5.0: Truncate jacobian integrals?" ,false,"Set to true for truncation of jacobian integrals below some threshhold.");
    Coarsening_List->set("C5.1: Truncation threshhold" ,0.0001,"Set cutoff threshold value for truncation of jacobian integrals.");

    Coarsening_List->set("C2.0: Type of Jacobian coarsening","none","Select a method for coarsening the mesh or jacobian",JacCoarsenTypeValidator);
    Coarsening_List->set("C2.1: Esize_Jac",0.25,"Set a fixed mesh spacing to be used in calcuation of Jacobian integrals");

    Coarsening_List->set("C3: Nzone",1,"Set number of coarsened zones in the problem (don't count the most refined zone).");

    Array<double> Rmin_zone_Array( (Coarsening_List->get<int>("C3: Nzone")),0.0);
    Coarsening_List->set("C4: Rmin for each zone", Rmin_zone_Array, str_zones);

    Coarsening_List->set("C6.0: 1D boundary zone?", false, "true if the problem has a 1D limit at the edge of the domain");
    Coarsening_List->set("C6.1: Dim_1D_bc", 0, "indicate the dimension (x=0,y=1,z=2) where 1D boundary zone should be applied",DimValidator);
    Coarsening_List->set("C6.2: X_1D_bc", 0.0, "set a distance away from domain boundaries where 1D boundary zone");

    PhysicsMethod_List->set("PM1: Attractions in A22 Block of matrix?", false, "Set to true to move attractions from A12 block of Schur matrix to A22 block of Shur Matrix.\n");
    PhysicsMethod_List->set("PM2: Physics Scaling?", "No Physics Scaling", "Set to true to turn on physics scaling.\n Physics scaling attempts to moderate the effect of terms like exp(mu) in the JDC polymer options",PhysicsScalingValidator);
    PhysicsMethod_List->set("PM3: Analytic Jacobian?", true, "Set to true for analytic Jacobian.\n Set to false for a physics based approximate Jacobian for JDC polymers");


    Array<double> ScaleFac_Array(Fluid_List->get<int>("F1_Ncomp"),1.0);
    PhysicsMethod_List->set("PM4: ScaleFac[icomp]", ScaleFac_Array, str_scalefac);


    NonlinearSolver_List->set("NLS1: Nonlinear Solver", "Newton Built-In", "Select nonlinear solver type you would like to use.", NLSolverTypeValidator);
    NonlinearSolver_List->set("NLS2: Max Nonlinear Iterations", 10, "Set the Maximum number of nonlinear iterations", Niter_Validator);
    NonlinearSolver_List->set("NLS3: Newton Tolerance Absolute", 0.00000001, "Set absolute convergence tolerance for nonlinear solver",NLtolValidator);
    NonlinearSolver_List->set("NLS4: Newton Tolerance Relative", 0.0001, "Set relative convergence tolerance for nonlinear solver",NLtolValidator);
    NonlinearSolver_List->set("NLS5: Minimum update fraction", 1.0, str_updatefrac,NLupdate_Validator);
    NonlinearSolver_List->set("NLS6: Picard Tolerance Absolute", 0.00000001, "Set absolute convergence tolerance for nonlinear solver",NLtolValidator);
    NonlinearSolver_List->set("NLS7: Picard Tolerance Relative", 0.0001, "Set relative convergence tolerance for nonlinear solver",NLtolValidator);

    LinearSolver_List->set("LS1: Linear Solver", "Schur Solver", "Select linear solver type you would like to use.",LSolverTypeValidator);
    LinearSolver_List->set("LS2: Max Linear Iterations", 100, "Set the Maximum number of nonlinear iterations", Niter_Validator);
    LinearSolver_List->set("LS3: Lin. Solver Tolerance", 0.001, "Set convergence tolerance for linear solver");
    LinearSolver_List->set("LS4: Linear Solver Approach", "GMRES", "Select linear solver approach you would like to use.",LSolverApproachValidator);
    LinearSolver_List->set("LS5: Matrix Scaling option", "none", "Select linear solver approach you would like to use.",LScalingValidator);
    LinearSolver_List->set("LS6: Preconditioner option", "none", "Select preconditioner approach you would like to use.",LPreconditionerValidator);
    LinearSolver_List->set("LS7: Number of Fill Levels for ILUT", 4.0, "Set the levels for ilut preconditioner");

  return;
}
/*****************************************************************************************************************************/
void dft_GUI_NumericalMethods_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                     Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                     Teuchos::RCP<Teuchos::ParameterList> Solver_List,
                     Teuchos::RCP<Teuchos::ParameterList> Coarsening_List,
                     Teuchos::RCP<Teuchos::ParameterList> LoadBalance_List,
                     Teuchos::RCP<Teuchos::ParameterList> PhysicsMethod_List,
                     Teuchos::RCP<Teuchos::ParameterList> NonlinearSolver_List,
                     Teuchos::RCP<Teuchos::ParameterList> LinearSolver_List)
{
  string str_coarse,str_zones,str_updatefrac,str_scalefac;
  int i,j;
    /**************************************/
    /* STRINGS */
    /**************************************/
     str_coarse="Select a method for coarsening of the problem.\n options: bulk_zone: Set to true to replace DFT Euler-Lagrange equation with rho(r)=rho_b some distance from the surfaces.\n \t Poisson-Bolzmann zone: Set to true to replace DFT Euler-Lagrange equation with rho(r)=rho_b some distance from the surfaces.\n\t Residual coarsening: Coarsening of residual equations and Jacobian integrals in zones away from the surfaces.\n\t Jacobian coarsening: Coarsening of jacobian integrals only in zones away from surfaces.";

    str_zones= "Set the minimum distance from surface in each zone. \n The distances should be arranged from nearest to furthest where\n the nearest zone to the surface is the most refined zone,\n and the furthest zone from the surfaces is the most coarse zone.  The first entry should be 0.0";

    str_updatefrac= "Set a minimum fraction for solution updates if using the built in Newton solver.\n This allows code to perform fractional updates if convergence difficulties are identifed. Set to 1.0 for full Newton steps. \n For Picard solver, this fraction is a mixing parameter for old and new solutions.  If using the mixed Picard-Newton Solver enter the minimum Newton update fraction.\n The Picard mixing parameter will be set to the newton update parameter divided by 100.";

    str_scalefac="Enter the desired scaling factor for segment component on each polymer chain.  \nNote that there will be zero entries wherever there are no segments of a given type on a particular chain";
    /**************************************/
    /* VALIDATORS */
    /**************************************/
 
 
    RCP<StringValidator> NLSolverTypeValidator = rcp(
           new StringValidator(tuple<std::string>("Newton Built-In","Newton NOX", 
						  "Picard Built-In","Picard NOX",
						  "Picard/Newton Built-In","Picard/Newton NOX")));

    RCP<StringValidator> PhysicsScalingValidator = rcp(
           new StringValidator(tuple<std::string>("No Physics Scaling","Automatic Calculation","Manual Input")));

    RCP<StringValidator> CoarsenTypeValidator = rcp(
           new StringValidator(tuple<std::string>("none","Bulk zone","Poisson-Boltzmann zone",
                                  "residual and jacobian coarsening","jacobian only coarsening")));

    RCP<EnhancedNumberValidator<int> > Niter_Validator = rcp(new EnhancedNumberValidator<int>());

    RCP<EnhancedNumberValidator<double> > NLtolValidator = rcp(new EnhancedNumberValidator<double>(1.e-10,1.e-2,1.e-6,8));

    RCP<EnhancedNumberValidator<double> > NLupdate_Validator = rcp(new EnhancedNumberValidator<double>(0.001,1.0,0.01));

    RCP<StringValidator> LoadBalValidator = rcp(
           new StringValidator(tuple<std::string>("Linear Matrix Balance","Geometric Recursive Bisection","Weighted Recursive Bisection")));

    RCP<StringValidator> LoadBalWeightValidator = rcp(
           new StringValidator(tuple<std::string>("Set automatically","Make all weights equal")));

    RCP<StringValidator> LSolverTypeValidator = rcp(
           new StringValidator(tuple<std::string>("Schur Solver","General/AztecOO solver")));

    RCP<StringValidator> LSolverApproachValidator = rcp(
           new StringValidator(tuple<std::string>("GMRES","cg","tfqmr","cgs","bicgstab")));

    RCP<StringValidator> LScalingValidator = rcp(
           new StringValidator(tuple<std::string>("none","row_sum","jacobi","symrow_sum")));

    RCP<StringValidator> LPreconditionerValidator = rcp(
           new StringValidator(tuple<std::string>("none","ilu","ilut","Jacobi","symmetric Gauss-Seidel","LSpoly3")));

    RCP<EnhancedNumberValidator<int> > DimValidator = rcp(new EnhancedNumberValidator<int>(0,2,1));

    RCP<StringValidator> JacCoarsenTypeValidator = rcp(
         new StringValidator(tuple<std::string>("none",
        "Jacobian Coarsening identical to Residual Coarsening",
	"Jacobian: factor of 2 more coarse than resid in most refined zone",
	"Jacobian: factor of 2 more coarse than resid in all but most coarse zone",
	"Jacobian: use most coarse zone to define entire matrix",
	"Jacobian: use 2nd most coarse zone in all but most coarse region",
	"Jacobian: set Esize_jac for all matrix calculations")));

   if (Load_Bal_Flag==LB_LINEAR) 
        LoadBalance_List->set("LB1: Load Balancing Approach", "Linear Matrix Balance", "Select a method for load balancing a parallel processing job", LoadBalValidator);
   else if (Load_Bal_Flag==LB_WEIGHTS) 
        LoadBalance_List->set("LB1: Load Balancing Approach", "Weighted Recursive Bisection", "Select a method for load balancing a parallel processing job", LoadBalValidator);
   else if (Load_Bal_Flag==LB_BOX) 
        LoadBalance_List->set("LB1: Load Balancing Approach", "Geometric Recursive Bisection", "Select a method for load balancing a parallel processing job", LoadBalValidator);

    LoadBalance_List->set("LB2: Load Balancing Weights", "Set automatically", "How to assign weights for weighted rcb method", LoadBalWeightValidator);

   if (Mesh_coarsening==FALSE)          Coarsening_List->set("C1: Coarsening Type", "none", str_coarse, CoarsenTypeValidator);
   else if (Mesh_coarsening==TRUE)      Coarsening_List->set("C1: Coarsening Type", "residual and jacobian coarsening", str_coarse, CoarsenTypeValidator);
   else if (Mesh_coarsening==PB_ZONE)   Coarsening_List->set("C1: Coarsening Type", "Poisson-Boltzmann zone", str_coarse, CoarsenTypeValidator);
   else if (Mesh_coarsening==BULK_ZONE) Coarsening_List->set("C1: Coarsening Type", "Bulk zone", str_coarse, CoarsenTypeValidator);


   if (Lcut_jac==TRUE){
      Coarsening_List->set("C5.0: Truncate jacobian integrals?" ,true,"Set to true for truncation of jacobian integrals below some threshhold.");
      Coarsening_List->set("C5.1: Truncation threshhold" ,Jac_threshold,"Set cutoff threshold value for truncation of jacobian integrals.");
   }

   if (Coarser_jac==JAC_RESID_ZONES_SAME) 
      Coarsening_List->set("C2.0: Type of Jacobian coarsening","Jacobian Coarsening identical to Residual Coarsening",
                                      "Select a method for coarsening the mesh or jacobian",JacCoarsenTypeValidator);
   else if (Coarser_jac==JAC_ZONE0_FAC2LESSTHANRESID)
      Coarsening_List->set("C2.0: Type of Jacobian coarsening","Jacobian: factor of 2 more coarse than resid in most refined zone",
                                      "Select a method for coarsening the mesh or jacobian",JacCoarsenTypeValidator);
   else if (Coarser_jac==JAC_ZONES_FAC2LESSTHANRESID)
      Coarsening_List->set("C2.0: Type of Jacobian coarsening","Jacobian: factor of 2 more coarse than resid in all but most coarse zone",
                                      "Select a method for coarsening the mesh or jacobian",JacCoarsenTypeValidator);
   else if (Coarser_jac==JAC_ZONES_ALLMOSTCOARSE)
      Coarsening_List->set("C2.0: Type of Jacobian coarsening","Jacobian: use most coarse zone to define entire matrix",
                                      "Select a method for coarsening the mesh or jacobian",JacCoarsenTypeValidator);
   else if (Coarser_jac==JAC_ZONES_SECONDMOSTCOARSE)
      Coarsening_List->set("C2.0: Type of Jacobian coarsening","Jacobian: use 2nd most coarse zone in all but most coarse region",
                                      "Select a method for coarsening the mesh or jacobian",JacCoarsenTypeValidator);
   else if (Coarser_jac==JAC_ZONES_SETFIXED_ESIZE)  
      Coarsening_List->set("C2.0: Type of Jacobian coarsening","Jacobian: set Esize_jac for all matrix calculations",
                                      "Select a method for coarsening the mesh or jacobian",JacCoarsenTypeValidator);

   Coarsening_List->set("C2.1: Esize_Jac",Jac_grid,"Set a fixed mesh spacing to be used in calcuation of Jacobian integrals");

   Coarsening_List->set("C3: Nzone",Nzone,"Set number of coarsened zones in the problem (don't count the most refined zone).");

   Array<double> Rmin_zone_Array(Nzone,0.0);
   if (Nzone > 1){
      for (i=0;i<Nzone;i++){
         if (i==0) Rmin_zone_Array[i]=0.0;
         else Rmin_zone_Array[i]=Rmax_zone[i-1];
      }
      Coarsening_List->set("C4: Rmin for each zone", Rmin_zone_Array, str_zones);
   }

   if (L1D_bc==TRUE){
     Coarsening_List->set("C6.0: 1D boundary zone?", true, "true if the problem has a 1D limit at the edge of the domain");
     Coarsening_List->set("C6.1: Dim_1D_bc", Grad_dim, "indicate the dimension (x=0,y=1,z=2) where 1D boundary zone should be applied",DimValidator);
     Coarsening_List->set("C6.2: X_1D_bc", X_1D_bc, "set a distance away from domain boundaries where 1D boundary zone will be applied");
   }

   if (ATTInA22Block==TRUE) 
       PhysicsMethod_List->set("PM1: Attractions in A22 Block of matrix?", true, "Set to true to move attractions from A12 block of Schur matrix to A22 block of Shur Matrix.\n");

   if (Physics_scaling==FALSE) 
       PhysicsMethod_List->set("PM2: Physics Scaling?", "No Physics Scaling", "Set to true to turn on physics scaling.\n Physics scaling attempts to moderate the effect of terms like exp(mu) in the JDC polymer options",PhysicsScalingValidator);
   else if (Physics_scaling==AUTOMATIC)
       PhysicsMethod_List->set("PM2: Physics Scaling?", "Automatic Calculation", "Set to true to turn on physics scaling.\n Physics scaling attempts to moderate the effect of terms like exp(mu) in the JDC polymer options",PhysicsScalingValidator);
   else if (Physics_scaling==MANUAL_INPUT) 
       PhysicsMethod_List->set("PM2: Physics Scaling?", "Manual Input", "Set to true to turn on physics scaling.\n Physics scaling attempts to moderate the effect of terms like exp(mu) in the JDC polymer options",PhysicsScalingValidator);
   if (Analyt_WJDC_Jac==FALSE) 
       PhysicsMethod_List->set("PM3: Analytic Jacobian?", false, "Set to true for analytic Jacobian.\n Set to false for a physics based approximate Jacobian for JDC polymers");

   Array<double> ScaleFac_Array(Ncomp,1.0);
   for (i=0; i<Npol_comp;i++) 
       for (j=0; j<Nblock[i];j++){
            ScaleFac_Array[SegType_per_block[i][j]]=Scale_fac_WJDC[i][SegType_per_block[i][j]];
   }
   PhysicsMethod_List->set("PM4: ScaleFac[icomp]", ScaleFac_Array, str_scalefac);

   if (NL_Solver==NEWTON_BUILT_IN)
      NonlinearSolver_List->set("NLS1: Nonlinear Solver", "Newton Built-In", "Select nonlinear solver type you would like to use.", NLSolverTypeValidator);
   else if (NL_Solver==NEWTON_NOX)
      NonlinearSolver_List->set("NLS1: Nonlinear Solver", "Newton NOX", "Select nonlinear solver type you would like to use.", NLSolverTypeValidator);
   else if (NL_Solver==PICARD_BUILT_IN)
      NonlinearSolver_List->set("NLS1: Nonlinear Solver", "Picard Built-In", "Select nonlinear solver type you would like to use.", NLSolverTypeValidator);
   else if (NL_Solver==PICARD_NOX)
      NonlinearSolver_List->set("NLS1: Nonlinear Solver", "Picard NOX", "Select nonlinear solver type you would like to use.", NLSolverTypeValidator);
   else if (NL_Solver==PICNEWTON_BUILT_IN)
      NonlinearSolver_List->set("NLS1: Nonlinear Solver", "Picard/Newton Built-In", "Select nonlinear solver type you would like to use.", NLSolverTypeValidator);
   else if (NL_Solver==PICNEWTON_NOX)
      NonlinearSolver_List->set("NLS1: Nonlinear Solver", "Picard/Newton NOX", "Select nonlinear solver type you would like to use.", NLSolverTypeValidator);

   NonlinearSolver_List->set("NLS2: Max Nonlinear Iterations", Max_NL_iter, "Set the Maximum number of nonlinear iterations", Niter_Validator);
   NonlinearSolver_List->set("NLS3: Newton Tolerance Absolute", NL_abs_tol, "Set absolute convergence tolerance for nonlinear solver",NLtolValidator);
   NonlinearSolver_List->set("NLS4: Newton Tolerance Relative", NL_rel_tol, "Set relative convergence tolerance for nonlinear solver",NLtolValidator);
   NonlinearSolver_List->set("NLS5: Minimum update fraction", NL_update_scalingParam, str_updatefrac,NLupdate_Validator);
   NonlinearSolver_List->set("NLS6: Picard Tolerance Absolute", NL_abs_tol_picard, "Set absolute convergence tolerance for nonlinear solver",NLtolValidator);
   NonlinearSolver_List->set("NLS7: Picard Tolerance Relative", NL_rel_tol_picard, "Set relative convergence tolerance for nonlinear solver",NLtolValidator);

   if (L_Schur==FALSE) LinearSolver_List->set("LS1: Linear Solver", "General/AztecOO solver", "Select linear solver type you would like to use.",LSolverTypeValidator);
   LinearSolver_List->set("LS2: Max Linear Iterations", Max_gmres_iter, "Set the Maximum number of linear iterations", Niter_Validator);
   LinearSolver_List->set("LS3: Lin. Solver Tolerance", Az_tolerance, "Set convergence tolerance for linear solver");

   if (Az_solver==0)
      LinearSolver_List->set("LS4: Linear Solver Approach", "GMRES", "Select linear solver approach you would like to use.",LSolverApproachValidator);
   else if (Az_solver==1)
      LinearSolver_List->set("LS4: Linear Solver Approach", "cg", "Select linear solver approach you would like to use.",LSolverApproachValidator);
   else if (Az_solver==2)
      LinearSolver_List->set("LS4: Linear Solver Approach", "tfqmr", "Select linear solver approach you would like to use.",LSolverApproachValidator);
   else if (Az_solver==3)
      LinearSolver_List->set("LS4: Linear Solver Approach", "cgs", "Select linear solver approach you would like to use.",LSolverApproachValidator);
   else if (Az_solver==4)
      LinearSolver_List->set("LS4: Linear Solver Approach", "bicgstab", "Select linear solver approach you would like to use.",LSolverApproachValidator);

   if (Az_scaling==-1)
       LinearSolver_List->set("LS5: Matrix Scaling option", "none", "Select linear solver approach you would like to use.",LScalingValidator);
   else if (Az_scaling==0)
       LinearSolver_List->set("LS5: Matrix Scaling option", "row_sum", "Select linear solver approach you would like to use.",LScalingValidator);
   else if (Az_scaling==1)
       LinearSolver_List->set("LS5: Matrix Scaling option", "jacobi", "Select linear solver approach you would like to use.",LScalingValidator);
   else if (Az_scaling==2)
       LinearSolver_List->set("LS5: Matrix Scaling option", "symrow_sum", "Select linear solver approach you would like to use.",LScalingValidator);

   if (Az_preconditioner==-1)
      LinearSolver_List->set("LS6: Preconditioner option", "none", "Select preconditioner approach you would like to use.",LPreconditionerValidator);
   else if (Az_preconditioner==0)
      LinearSolver_List->set("LS6: Preconditioner option", "ilu", "Select preconditioner approach you would like to use.",LPreconditionerValidator);
   else if (Az_preconditioner==1)
      LinearSolver_List->set("LS6: Preconditioner option", "Jacobi", "Select preconditioner approach you would like to use.",LPreconditionerValidator);
   else if (Az_preconditioner==2)
      LinearSolver_List->set("LS6: Preconditioner option", "symmetric Gauss-Seidel", "Select preconditioner approach you would like to use.",LPreconditionerValidator);
   else if (Az_preconditioner==3)
      LinearSolver_List->set("LS6: Preconditioner option", "LSpoly3", "Select preconditioner approach you would like to use.",LPreconditionerValidator);
   else if (Az_preconditioner==4)
      LinearSolver_List->set("LS6: Preconditioner option", "ilut", "Select preconditioner approach you would like to use.",LPreconditionerValidator);

   LinearSolver_List->set("LS7: Number of Fill Levels for ILUT", Az_ilut_fill_param, "Set the levels for ilut preconditioner");
   return;
}
/***************************************************************************************************************************************************/
void dft_GUI_NumericalMethods_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Polymer_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Solver_List,
                         Teuchos::RCP<Teuchos::ParameterList> Coarsening_List,
                         Teuchos::RCP<Teuchos::ParameterList> LoadBalance_List,
                         Teuchos::RCP<Teuchos::ParameterList> PhysicsMethod_List,
                         Teuchos::RCP<Teuchos::ParameterList> NonlinearSolver_List,
                         Teuchos::RCP<Teuchos::ParameterList> LinearSolver_List)
{

          /* dependencies for mesh/jacobian coarsening options */
     Dependency::ParameterEntryList Coarsening_Deps;
     Coarsening_Deps.insert(Coarsening_List->getEntryRCP("C2.0: Type of Jacobian coarsening"));
     Coarsening_Deps.insert(Coarsening_List->getEntryRCP("C3: Nzone"));
     Coarsening_Deps.insert(Coarsening_List->getEntryRCP("C4: Rmin for each zone"));
     RCP<StringVisualDependency> TypeCoarse_Dep = rcp(
           new StringVisualDependency(Coarsening_List->getEntryRCP("C1: Coarsening Type"),Coarsening_Deps,
               tuple<std::string>("residual and jacobian coarsening","jacobian only coarsening","Bulk zone","Poisson-Boltzmann zone")));

     Dependency::ParameterEntryList OneD_BoundaryZoneDeps;
     OneD_BoundaryZoneDeps.insert(Coarsening_List->getEntryRCP("C6.1: Dim_1D_bc"));
     OneD_BoundaryZoneDeps.insert(Coarsening_List->getEntryRCP("C6.2: X_1D_bc"));
     RCP<BoolVisualDependency> OneD_Boundary_Dep = rcp(
        new BoolVisualDependency(Coarsening_List->getEntryRCP("C6.0: 1D boundary zone?"),OneD_BoundaryZoneDeps,true));

     RCP<BoolVisualDependency> TruncThresh_Dep = rcp(
        new BoolVisualDependency(Coarsening_List->getEntryRCP("C5.0: Truncate jacobian integrals?"),Coarsening_List->getEntryRCP("C5.1: Truncation threshhold"),true));

      RCP<NumberArrayLengthDependency<int,double> > RminZoneLength_Dep = rcp(
          new NumberArrayLengthDependency<int,double>(Coarsening_List->getEntryRCP("C3: Nzone"), Coarsening_List->getEntryRCP("C4: Rmin for each zone")));

      RCP<StringVisualDependency> JacGrid_Dep = rcp(
           new StringVisualDependency(Coarsening_List->getEntryRCP("C2.0: Type of Jacobian coarsening"),
                                      Coarsening_List->getEntryRCP("C2.1: Esize_Jac"),"Jacobian: set Esize_jac for all matrix calculations"));

          /* dependencies for load balance options */

      RCP<StringVisualDependency> LoadBal_Dep = rcp(
           new StringVisualDependency(LoadBalance_List->getEntryRCP("LB1: Load Balancing Approach"),
                                      LoadBalance_List->getEntryRCP("LB2: Load Balancing Weights"),"Weighted Recursive Bisection"));

          /* dependencies for nonlinear solver options */
      Dependency::ParameterEntryList PicardSolverDeps;
      PicardSolverDeps.insert(NonlinearSolver_List->getEntryRCP("NLS6: Picard Tolerance Absolute"));
      PicardSolverDeps.insert(NonlinearSolver_List->getEntryRCP("NLS7: Picard Tolerance Relative"));
      RCP<StringVisualDependency> PicardSolver_Dep = rcp(
       new StringVisualDependency(NonlinearSolver_List->getEntryRCP("NLS1: Nonlinear Solver"),PicardSolverDeps,
           tuple<std::string>("Picard Built-In","Picard NOX", "Picard/Newton Built-In","Picard/Newton NOX")));

      Dependency::ParameterEntryList NewtonSolverDeps;
      NewtonSolverDeps.insert(NonlinearSolver_List->getEntryRCP("NLS3: Newton Tolerance Absolute"));
      NewtonSolverDeps.insert(NonlinearSolver_List->getEntryRCP("NLS4: Newton Tolerance Relative"));
      RCP<StringVisualDependency> NewtonSolver_Dep = rcp(
       new StringVisualDependency(NonlinearSolver_List->getEntryRCP("NLS1: Nonlinear Solver"),NewtonSolverDeps,
           tuple<std::string>("Newton Built-In","Newton NOX", "Picard/Newton Built-In","Picard/Newton NOX")));

          /* dependencies for physics methods solver options */
      RCP<StringVisualDependency> PhysATTA22_Dep = rcp(
           new StringVisualDependency(Functional_List->getEntryRCP("F2_PAIRPOTcore_Functional"),
                PhysicsMethod_List->getEntryRCP("PM1: Attractions in A22 Block of matrix?"), "No Mean Field Functional",false));

      RCP<StringCondition> PolyFuncCon = rcp(
           new StringCondition(Functional_List->getEntryRCP("F4_POLYMER_Functional"),
                              tuple<std::string>("Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)",
                              "Polymer_JDC_iSAFT(comp)")));

      RCP<StringCondition> ScaleFacManualCon = rcp(
           new StringCondition(PhysicsMethod_List->getEntryRCP("PM2: Physics Scaling?"),"Manual Input"));

      Condition::ConstConditionList ScaleFac_conList=tuple<RCP<const Condition> >(PolyFuncCon,ScaleFacManualCon);
      RCP<AndCondition> ScaleFac_Con = rcp(new AndCondition(ScaleFac_conList));
        RCP<ConditionVisualDependency> ScaleFac_Dep = rcp(
            new ConditionVisualDependency(ScaleFac_Con, PhysicsMethod_List->getEntryRCP("PM4: ScaleFac[icomp]")));

      Dependency::ParameterEntryList ScaleFacArrayLength_Deps;
      ScaleFacArrayLength_Deps.insert(PhysicsMethod_List->getEntryRCP("PM4: ScaleFac[icomp]"));
      RCP<NumberArrayLengthDependency<int,double> > ScaleFacLength_Dep = rcp(
           new NumberArrayLengthDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"), ScaleFacArrayLength_Deps));

/*      Dependency::ParameterEntryList ScaleFac_Deps;
      ScaleFac_Deps.insert(PhysicsMethod_List->getEntryRCP("PM4: ScaleFac[ipol][icomp]"));
          RCP<TwoDRowDependency<int,double> > ScaleFacRows_Dep = rcp(
           new TwoDRowDependency<int,double>(Polymer_List->getEntryRCP("P1: Npoly_comp"), ScaleFac_Deps));

          RCP<TwoDColDependency<int,double> > ScaleFacCol_Dep = rcp(
           new TwoDColDependency<int,double>(Fluid_List->getEntryRCP("F1_Ncomp"), ScaleFac_Deps));*/


      RCP<StringVisualDependency> PhysAnalytJac_Dep = rcp(
           new StringVisualDependency(Functional_List->getEntryRCP("F4_POLYMER_Functional"),PhysicsMethod_List->getEntryRCP("PM2: Physics Scaling?"),
               tuple<std::string>("Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)","Polymer_JDC_iSAFT(comp)")));

      RCP<StringVisualDependency> PhysScale_Dep = rcp(
           new StringVisualDependency(Functional_List->getEntryRCP("F4_POLYMER_Functional"),PhysicsMethod_List->getEntryRCP("PM3: Analytic Jacobian?"),
               tuple<std::string>("Polymer_JDC_iSAFT(seg)","Polymer_JDC_iSAFT(segRho compField)","Polymer_JDC_iSAFT(comp)")));

          /* dependencies for linear solver options */

      Dependency::ParameterEntryList GenSolverDeps;
      GenSolverDeps.insert(LinearSolver_List->getEntryRCP("LS4: Linear Solver Approach"));
      GenSolverDeps.insert(LinearSolver_List->getEntryRCP("LS5: Matrix Scaling option"));
      GenSolverDeps.insert(LinearSolver_List->getEntryRCP("LS6: Preconditioner option"));
      RCP<StringVisualDependency> GenSolver_Dep = rcp(
       new StringVisualDependency(LinearSolver_List->getEntryRCP("LS1: Linear Solver"),GenSolverDeps,"General/AztecOO solver"));

      RCP<StringVisualDependency> LevelILUT_Dep = rcp(
           new StringVisualDependency( LinearSolver_List->getEntryRCP("LS6: Preconditioner option"),
                                       LinearSolver_List->getEntryRCP("LS7: Number of Fill Levels for ILUT"), "ilut"));

    /*****************************************/
    /* add the dependencies for this section.*/
    /*****************************************/
      depSheet_Tramonto->addDependency(TruncThresh_Dep);
      depSheet_Tramonto->addDependency(OneD_Boundary_Dep);
      depSheet_Tramonto->addDependency(TypeCoarse_Dep);
      depSheet_Tramonto->addDependency(RminZoneLength_Dep);
      depSheet_Tramonto->addDependency(JacGrid_Dep);
      depSheet_Tramonto->addDependency(LoadBal_Dep);
      depSheet_Tramonto->addDependency(PicardSolver_Dep);
      depSheet_Tramonto->addDependency(NewtonSolver_Dep);
      depSheet_Tramonto->addDependency(PhysATTA22_Dep);
      depSheet_Tramonto->addDependency(PhysAnalytJac_Dep);
      depSheet_Tramonto->addDependency(PhysScale_Dep);
      depSheet_Tramonto->addDependency(GenSolver_Dep);
      depSheet_Tramonto->addDependency(LevelILUT_Dep);
      depSheet_Tramonto->addDependency(ScaleFac_Dep);
      depSheet_Tramonto->addDependency(ScaleFacLength_Dep);
      /*depSheet_Tramonto->addDependency(ScaleFacRows_Dep);
      depSheet_Tramonto->addDependency(ScaleFacCol_Dep);*/

    return;
}

