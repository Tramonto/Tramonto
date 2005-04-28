/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
 
/*
 *  FILE: dft_newton.c
 *
 *  This file solves the dft problem using Newton's method. The fill
 *  of the residual and matrix are in file dft_fill.c, and the linear
 *  solve is done using calls to the Aztec library.
 */
#include "mpi.h" 


#ifdef HAVE_DFT_SCHUR_SOLVER
#include "dft_schur_solver.h"
#endif

#include "dft_globals_const.h" 
#include "rf_allo.h"

/* NEWTON SOLVER using dft_SolverManager */
int solve_problem(double **x_internal_ptr, double **x2_internal_ptr)
/*
 * x and x2 in dfrt_main have been renamed  "internal" to represent 1d arrays
 * x2 only relevent for Lbinodal calculations, so can mostly be ignored
 */

{
  int iter;
  double **x, **x2, **xLocal;

  /* Construct dft_SolverManager with information on number of unknowns*/
  Solver_manager = dft_solvermanager_create(Nunk_per_node, Unk2Phys, NULL, NULL, MPI_COMM_WORLD);

  /* Give Nodal Row and Column maps */
  (void) dft_solvermanager_setnodalrowmap(Solver_manager, Nnodes_per_proc, L2G_node);
  (void) dft_solvermanager_setnodalcolmap(Solver_manager, Nnodes_box     , B2G_node);

  /* SolverManager can now set up its own numbering scheme, set up unknown-based Maps */
  (void) dft_solvermanager_finalizeblockstructure(Solver_manager);

  /* Create 2d arrays for solution vector and reconcile ghost nodes using import */
  xLocal = (double **) array_alloc(2, Nunk_per_node, Nnodes_per_proc, sizeof(double));
  x      = (double **) array_alloc(2, Nunk_per_node, Nnodes_box     , sizeof(double));
  internal2local(*x_internal_ptr, xLocal);
  (void) dft_solvermanager_importR2C(Solver_manager, xLocal, x);

  /* Same for second solution vector when Lbinodal is true */
  if (Lbinodal) {
    x2   = (double **) array_alloc(2, Nunk_per_node, Nnodes_box     , sizeof(double));
    internal2local(*x2_internal_ptr, xLocal);
    (void) dft_solvermanager_importR2C(Solver_manager, xLocal, x2);
  }

#ifdef LOCA_REIMPLEMENTED
  if (Loca.method != -1)
    iter = solve_continuation(x, x2, Sten_Type[POLYMER_CR], aux_info);
  else
#endif
     iter = newton_solver(x, NULL);

  /* Translate solution back into internal (1d) coordinates */
  
  box2local(x, xLocal);
  loca2internal(xLocal, *x_internal_ptr);

  if (Lbinodal) {
    box2local(x2, xLocal);
    loca2internal(xLocal, *x2_internal_ptr);
  }

  /* Call the destructor for the dft_SolverManager */
  safe_free((void **) &xLocal);
  dft_solvermanager_destruct(Solver_manager);

  return(iter);
}

/****************************************************************************/
/****************************************************************************/
void internal2local(double* xInternal, double** xLocal) {
  /* Routine to translate 1d array of owned unknows to 2d array */
  int iunk, inode;

  for (iunk=0; iunk<Nunk_per_node; iunk++)
    for (inode=0; inode<Nnodes_per_proc; inode++)
      xLocal[iunk][inode] = xInternal[loc_find(iunk,inode,LOCAL)];
}

/****************************************************************************/
/****************************************************************************/
void local2internal(double** xLocal, double* xInternal) {
  /* Routine to translate 1d array of owned unknows to 2d array */
  int iunk, inode;

  for (iunk=0; iunk<Nunk_per_node; iunk++)
    for (inode=0; inode<Nnodes_per_proc; inode++)
      xInternal[loc_find(iunk,inode,LOCAL)] = xLocal[iunk][inode];
}

/****************************************************************************/
/****************************************************************************/
void box2local(double** xBox, double** xLocal) {
  /* Routine to restrict 2d box array to 2d node array */
  /* This is the opposite of importR2C, and strips ghost nodes */
  int iunk, inode, ibox;

  for (ibox=0; ibox<Nnodes_box; ibox++) {
    inode = B2L_node[ibox];
    if (inode!=-1) {
      for (iunk=0; iunk<Nunk_per_node; iunk++)
        xLocal[iunk][inode] = xBox[iunk][ibox];
    }
  }
}

/****************************************************************************/
/****************************************************************************/
int newton_solver(double** x, void* con_ptr) {

  int iter=0;
  int converged=FALSE, converged2=TRUE;
  double** delta_x;
  delta_x = (double **) array_alloc(2, Nunk_per_node, Nnodes_box, sizeof(double));

  do {
    iter++;

    (void) dft_solvermanager_initializeproblemvalues(Solver_manager);

    /* Call Matrix and Residual Fill routine, resid_only_flag=FALSE)*/
    fill_resid_and_matrix_control(x, FALSE);

    (void) dft_solvermanager_finalizeproblemvalues(Solver_manager);
    (void) dft_solvermanager_setupsolver(Solver_manager);
    (void) dft_solvermanager_solve(Solver_manager);
    
    /* I am assuming getLhs returns box coordinates (e.g. Column Map)!! */
    (void) dft_solvermanager_getlhs(Solver_manager, delta_x);
    
    /* Do: x += delta_x, and check for convergence */
    converged = update_solution(x, delta_x, iter);

  } while (iter < Max_Newton_iter && (!converged || !converged2));

  if (!converged || !converged2) {
    printf("\tNewton Solver: Failed to converge in %d iterations\n",iter);
    iter = -iter;
  }
  else
    printf("\tNewton Solver: Successful convergence in %d iterations\n",iter);

  safe_free((void *) &delta_x);
  return iter;
}

/****************************************************************************/
/****************************************************************************/
int update_solution(double** x, double** delta_x, int iter) {
/* Routine to update solution vector x using delta_x and 
 * to test for convergence of Newton's method.
 * 
 * Note: Modifications to Newton's method, including damping, have not yet
 *       been translated form previous update_solution.
 *       iter  value may be used in some damping methods
 */

  int iunk, ibox, inode;
  double updateNorm=0.0, temp;
  char *yo = "newupdate solution";

  for (ibox=0; ibox<Nnodes_box; ibox++) {

    /* Increment updateNorm only for owned nodes (inode=-1 for ghosts) */
    inode = B2L_node[ibox];
    if (inode != -1) {
      for (iunk=0; iunk<Nunk_per_node; iunk++) {
        temp = delta_x[iunk][ibox]/(Newton_rel_tol*x[iunk][ibox] + Newton_abs_tol);
        updateNorm +=  temp*temp;
      }
    }

    /* Update all solution componenets */
    for (iunk=0; iunk<Nunk_per_node; iunk++)
      x[iunk][ibox] += delta_x[iunk][ibox];
  }

  /* Aztec call can be replaced by MPI? */
  updateNorm = sqrt(gsum_double(updateNorm));

  if (Proc==0 && Iwrite != NO_SCREEN)
    printf("\t\t%s: Weighted norm of update vector =  %g\n", yo, updateNorm);

  if (updateNorm > 1.0) return(FALSE);
  else                  return(TRUE);

}
