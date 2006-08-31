/*
//@HEADER
// ********************************************************************
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER
*/

/*
Pseudo code for Tramonto interacting with modified dft_Solvermanager
Prototype code for [unknown#][node] data structures in Tramonto

NOTE: The term "Unknown" now means exclusively the undiscretized
quantity (e.g. density, rho_bar_0) described by the physics being
solved and not the discretized number of unknowns. This term did 
have two meanings, but in the data structure I am proposing here
there is no need to name the discretized degrees-of-freedom.
*/

///Tramonto Unknown Information:
int numUnknowns; //sum of previous array over all blocks
int Phys[numUnknowns]; //Identifier of Physics type, e.g DENSITY, POISSON, RHO_BARS,...

// Identify the unknowns being solved (can be switched to physics names, 
// instead of unknown names....
numUnknowns=7;

Phys[0]=PHO_BAR_0;
Phys[1]=PHO_BAR_1;
Phys[2]=PHO_BAR_2;
Phys[3]=PHO_BAR_3;
Phys[4]=PHO_BAR_V1;
Phys[5]=PHO_BAR_V2;
Phys[6]=DENSITY;
Phys[7]=POISSON;

//Tramonto Mesh Variables:
int numOwnedNodes, numBoxNodes;
  // Maps: Owned-to-global, box-to-global, local-to-box
int OwnedNodes[numOwnedNodes], BoxNodes[numBoxNodes], L2B[numOwnedNodes];

//Construct dft_Solvermanager: 
//Physics names passed in so that block decisions can be made in this class
dft_solvermanager_create(NumUnknowns, Phys, MPI_COMM_WORLD);

//Added "Nodal" to the name to highlight the change from previous version
dft_solvermanager_setNodalrowmap(numOwnedNodes, OwnedNodes);
dft_solvermanager_setNodalcolumnmap(numBoxNodes, BoxNodes);

// The following routine creates the real Maps for the matrix blocks
// but Tramonto doesn't need to know them. The choice of node ordering
// or physics Unknown numbering can be made internally here, as can
// the choice to form 1 or many epetra matrices, Crs vs VBR, etc.
dft_solvermanager_finalizeblockstructure();


// Mike: No need for Graph methods in dft_SolutionManager, we want
// to fill the matrix directly wiithout a graph-building stage.



// Pseudo code for matrix and residual fill:
dft_solvermanager_initializeproblemvalues();

for (inode=0; inode<numOwnedNodes; inode++) {
  // Residual vector loaded with row inode, solution vector
  // indexed with row jbox, diagonal entries in matrix
  // loaded withnode numberings  [inode][ibox]
  ibox = L2B(inode);

  // Loop over unknown types at each node. This loop could be changed
  // so that chunks of physics get loaded together, like all rho bars,
  for (iunk=0; iunk<numUnknowns; iunk++) {
    if (Phys[iunk]==POISSON) {
      for (junk=0; junk<numUnknowns; junk++) {
        if (Phys[junk]==POISSON) {
	  // Now doing residual cand jacobian contributions
	  // for Poisson EQ w.r.t Poisson Unknowns.

          // residual vector f aka Rhs is accessed with local nodes and iunk
          // solution vector x aka Lhs accessed with box nodes and junk
 
          if (FILL_RESIDUAL) {
            // Residual fill for 1D LaPlace operator f[iunk][inode]
            f = x[junk][ibox-1] -2*x[junk][ibox] + x[junk][ibox+1];
 
            //Row number is now expanded to 2 integers [unk#][node]
	    //I think we want this a sum-into method, not a replace!
            dft_solutionmanager_insertRhsValue(iunk, inode, -f);
          }
 
          if (FILL_JACOBIAN) {
            //Jacobian fill for same: this loads three entries in 1 call
            //nodeIndices is an integer temp array, Values a double temp array
            numEntries=3;
            nodeIndicies[0]=ibox-1; nodeIndicies[1]=ibox; nodeIndicies[2]=ibox+1; 
            Values[0]=1.0; Values[1]=-2.0; Values[2]=1.0; 
 
            // Row and Column are now each expressed with 2 integers: node,unknown
            // so an extra 2 arguments added to this routine: rowUnkNumber and colUnkNumber
            // I think the colUnkNumber can be an integer and not a more general array of 
            // integers.
            dft_solutionmanager_insertMatrixValues(iunk, inode, junk, numEntries, Values, nodeIndices);
          }
        } // Loop over junk==POISSON

	else if (Phys[junk]==DENSITY) {
	  // For the same equation, add contributions from density unknown
          if (FILL_RESIDUAL) {
            f = charge*factor*x[junk][ibox];
            dft_solutionmanager_insertRhsValue(iunk, inode, -f);
	  }
          if (FILL_JACOBIAN) {
            // Proposed single matrix entry interface -- one less argument, and value and indices are not pointers
            dft_solutionmanager_insertOneMatrixValue(iunk, inode, junk, factor*charge, ibox);
	  }
        }
      } //Loop over junk
    }//end IF POISSON EQ
    
    //All the above was for Phys[iunk]==POISSON. Now check for other physics

    // Some code fragments for rho_bar_0 fill w.r.t density unknowns
 
    if (Phys[iunk]==RHO_BAR_0) {
      for (junk=0; junk<numUnknowns; junk++) {
        if (Phys[junk]==DENSITY) {

          // Uses RB0 structure holding Stencil info for this integral
	  f=0;
          for (isten=0; isten<RB0.stenPoints; isten++) {
            jbox = RB0.StencilOffset(isten, ibox);
            weight = RB0.StencilWeight(isten);
            f -= weight*x[junk][jbox];
     
            // Add 1 entry into Jacobian
            dft_solutionmanager_insertOneMatrixValue(iunk, inode, junk, -weight, jbox);
          }
          // LOAD -Residual contribution from this node, unknown
          dft_solutionmanager_insertRhsValue(iunk, inode, -f);
        } // End of contributions of this equation due to density unknowns
      } //END of loop over unknowns 
    } //end of filling Rho_BAR_0 EQuation
  } //End of unknown iunk loop
}//end inode loop

// Mike: It would be great to have an internal flag e.g. "MatrixIsFinalized" 
// in dft_SolutionManager that allows us to fill the matrix graph/values
// from scratch the first fill, but then keeps the graph around after that.
// The first time this finalizeproblemvalues method gets called, OptimizeStorage
// gets called and then this new flag gets switched to true. I am not sure
// whether the  insertMatrixValues  routine will be different based on
// this flag (Insert for first fill vs SumInto for subsequent). A resetMatrix
// method can be added in case we want to change the sparsity pattern during 
// a run.
dft_solvermanager_finalizeproblemvalues();



//Solve/Newton
dft_solvermanager_solve();

// Newton update: Lhs is 2D array
dft_solvermanager_getLhs(delta_x);
for (i=0; i<numUnknowns; i++)
  for (j=0; j<numBoxNodes; j++)
    x[i][j] += delta_x[i][j];

printf("Total Time = %g (sec)\n", old_time*0.1);
