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
 *  FILE: dft_output.c
 *
 *  This file contains routines that post-process the density profiles.
 *
 */

#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

#define BFD 0
#define FFD 1
#define CFD 2

void sum_rho_wall (double *, double **);
void force_elec (double *, double **);
void find_pot_derivs (double *, double *);
double sum_rho_midplane (double *);
double calc_local_pressure (double *);
double calc_deriv (int, int, int, int *, double *, int);
void find_offset (int, int, int *);
void integrate_rho_vdash (double *, double **);

/**************************************************************************** 
 * calc_force: This routine calculates the contribution to the solvation 
               force in all directions on every wall for this processor */

void
calc_force (FILE * fp, double *x, double fac_area)
{
  double **p_tilde, **f_elec, area;
  double **p_tilde_L, **f_elec_L;
  double p_tilde_iwall_idim, f_elec_iwall_idim, force;
  int idim, iwall, i;

  p_tilde = (double **) array_alloc (2, Nwall, Ndim, sizeof (double));
  if (Ipot_wf_c == 1)
    f_elec = (double **) array_alloc (2, Nwall, Ndim, sizeof (double));

  if (!Lhard_surf && Restart != 4)
    integrate_rho_vdash (x, p_tilde);
  else
    {
      sum_rho_wall (x, p_tilde);
/*       rho_sum_mid = sum_rho_midplane(x);
         rho_sum_tot = AZ_gsum_double(rho_sum_mid,Aztec.proc_config);*/
    }

  /* calculate the electrostatic contribution to the force */
  if (Ipot_wf_c == 1)
    force_elec (x, f_elec);

  /* deal with linked surfaces */
  p_tilde_L = (double **) array_alloc (2, Nlink, Ndim, sizeof (double));
  if (Ipot_wf_c == 1)
    f_elec_L = (double **) array_alloc (2, Nlink, Ndim, sizeof (double));

  for (idim = 0; idim < Ndim; idim++)
    {
      for (i = 0; i < Nlink; i++)
	{
	  p_tilde_L[i][idim] = 0.0;
	  if (Ipot_wf_c == 1)
	    f_elec_L[i][idim] = 0.0;
	  for (iwall = 0; iwall < Nwall; iwall++)
	    {
	      if (Link[iwall] == i)
		{
		  p_tilde_L[i][idim] += p_tilde[iwall][idim];
		  if (Ipot_wf_c == 1)
		    f_elec_L[i][idim] += f_elec[iwall][idim];
		}
	    }
	}
    }

/* now collate all the contributions from each processor and then
   divide the sum total by the surface area of iwall */

  /* get sum of p_tilde[iwall][idim] and if you're proc #0 print it !! */
  if (Proc == 0 && Iwrite != NO_SCREEN)
    printf ("\n----------------------------------------------------------\n");

  for (i = 0; i < Nlink; i++)
    {
      for (idim = 0; idim < Ndim; idim++)
	{

	  f_elec_iwall_idim = 0.0;
	  if (Ipot_wf_c == 1)
	    f_elec_iwall_idim = AZ_gsum_double (f_elec_L[i][idim],
						Aztec.proc_config) *
	      Temp_elec / (4.0 * PI);

	  p_tilde_iwall_idim =
	    AZ_gsum_double (p_tilde_L[i][idim], Aztec.proc_config);

	  force = p_tilde_iwall_idim + f_elec_iwall_idim;

	  if (Proc == 0)
	    {

	      if (Lper_area)
		{
		  area = 0.0;
		  if (Nlink == Nwall)
		    area = S_area_tot[Nlists_HW - 1][0];
		  else
		    for (iwall = 0; iwall < Nwall; iwall++)
		      {
			if (Link[iwall] == 0)
			  area += S_area_tot[Nlists_HW - 1][iwall];
		      }
		  if (area > 0.0)
		    force /= area;
		  else
		    force *= fac_area;
		}
	      else
		{		/* not normalizing by area */
		  area = 1.0;
		  force *= fac_area;
		}

	      if (Iwrite != NO_SCREEN)
		{
		  printf ("iwall: %d \t idim: %d \n", i, idim);
		  printf ("fac: %9.6f  area: %9.6f\n", fac_area, area);
		  printf ("\t\t p_tilde[][]: %9.6f\n", p_tilde_iwall_idim);
		  printf ("\t\t f_elec [][]: %9.6f\n", f_elec_iwall_idim);
		  printf ("\t\t total force: %9.6f\n", force);
		}

	      if (i == 0 && idim == Orientation[WallType[0]])
		fprintf (fp, "%12.9f  ", force);
/*                   fprintf(fp,"%12.9f  %12.9f  %12.9f  ",
                          force,
                          p_tilde_iwall_idim/area,
                          f_elec_iwall_idim/area);*/
	    }
	}
    }


  if (Proc == 0 && Iwrite != NO_SCREEN)
    {
      printf ("----------------------------------------------------------\n");
    }
  safe_free ((void *) &p_tilde);
  safe_free ((void *) &p_tilde_L);
  if (Ipot_wf_c == 1)
    {
      safe_free ((void *) &f_elec);
      safe_free ((void *) &f_elec_L);
    }
  return;
}

/***************************************************************************
 * sum_rho_wall: For the hard wall case, sum the densities at the wall to  *
 *              get the force on the wall.                                 */

void
sum_rho_wall (double *x, double **Sum_rho)
{
  int loc_i, loc_inode, iwall, idim, icomp, ilist,
    iel_w, inode, surf_norm, inode_box, jwall;
  double prefac, nodepos[3];

  for (iwall = 0; iwall < Nwall; iwall++)
    {
      for (idim = 0; idim < Ndim; idim++)
	{
	  Sum_rho[iwall][idim] = 0.0;
	}
    }

  for (icomp = 0; icomp < Ncomp; icomp++)
    {
      if (Nlists_HW == 1 || Nlists_HW == 2)
	ilist = 0;
      else
	ilist = icomp;

      for (loc_inode = 0; loc_inode < Nnodes_per_proc; loc_inode++)
	{
	  inode = Aztec.update[Nunk_per_node * loc_inode] / Nunk_per_node;
	  inode_box = L2B_node[loc_inode];
	  iwall = Nodes_2_boundary_wall[ilist][inode_box];

	  if (iwall != -1)
	    {

	      node_to_position (inode, nodepos);

	      if (Sten_Type[POLYMER_CR])
		loc_i =
		  Aztec.update_index[Ncomp + icomp +
				     Nunk_per_node * loc_inode];
	      else
		loc_i = Aztec.update_index[icomp + Nunk_per_node * loc_inode];

	      for (iel_w = 0; iel_w < Nelems_S[ilist][loc_inode]; iel_w++)
		{

		  jwall = Surf_elem_to_wall[ilist][loc_inode][iel_w];
		  surf_norm = Surf_normal[ilist][loc_inode][iel_w];
		  idim = abs (surf_norm) - 1;
		  prefac = (double) (surf_norm / abs (surf_norm)) *
		    Area_surf_el[idim] / ((double) Nnodes_per_el_S);

		  /* taking mean between force limits from above (high rho)
		     and below (rho=0) require adding 1/2 the contributions
		     from above.  Remainder of profile is unaffected by 
		     these points so they have a weight of 1 */

		  if (iwall == -2)
		    Sum_rho[jwall][idim] += 0.5 * prefac * x[loc_i];
		  else
		    Sum_rho[jwall][idim] += prefac * x[loc_i];

		}		/* end of surface element loop */
	    }			/* end of test for boundary node */
	}			/* end of loop over nodes on this processor */
    }				/* end of icomp loop */
  return;
}

/***************************************************************************
 * force_elec: For charged surfaces, find the electrostatic contribution   *
 *             to the force on each wall.                                  */

void
force_elec (double *x, double **Sum_dphi_dx)
{
  int loc_i, loc_inode, iwall, idim, jdim,
    iel_w, inode, surf_norm, ijk_box[3],
    el_type, offset[3], loc_j, reflect_flag[3],
    match, test, inode_box, jnode_box, blocked, jwall;

  double prefac, nodepos[3],
    deriv_x[12][3], dot_prod[12], integrand, store1, store2,
    store1_tot, store2_tot;

  store1 = 0.0;
  store2 = 0.0;

  for (idim = 0; idim < Ndim; idim++)
    {
      for (iwall = 0; iwall < Nwall; iwall++)
	{
	  Sum_dphi_dx[iwall][idim] = 0.0;
	}
    }

  for (loc_inode = 0; loc_inode < Nnodes_per_proc; loc_inode++)
    {
      inode_box = L2B_node[loc_inode];
      inode = Aztec.update[Nunk_per_node * loc_inode] / Nunk_per_node;
      /*iwall = Nodes_2_boundary_wall[Nlists_HW-1][inode_box]; */
      iwall = Nodes_2_boundary_wall[0][inode_box];

      if (iwall != -1)
	{

	  node_box_to_ijk_box (inode_box, ijk_box);

	  node_to_position (inode, nodepos);

/*          for (iel_w=0; iel_w<Nelems_S[Nlists_HW-1][loc_inode]; iel_w++){
              surf_norm = Surf_normal[Nlists_HW-1][loc_inode][iel_w];*/
	  for (iel_w = 0; iel_w < Nelems_S[0][loc_inode]; iel_w++)
	    {
	      surf_norm = Surf_normal[0][loc_inode][iel_w];
	      if (iel_w == 0)
		{
		  test = surf_norm;
		  match = TRUE;
		}
	      else if (surf_norm != test)
		match = FALSE;
	    }

	  /* for (iel_w=0; iel_w<Nelems_S[Nlists_HW-1][loc_inode]; iel_w++){ */
	  for (iel_w = 0; iel_w < Nelems_S[0][loc_inode]; iel_w++)
	    {
	      dot_prod[iel_w] = 0.0;
	      for (idim = 0; idim < Ndim; idim++)
		deriv_x[iel_w][idim] = 0.0;
	    }

	  if (Matrix_fill_flag < 3)
	    loc_i = Aztec.update_index[Ncomp + Nunk_per_node * loc_inode];
	  else
	    loc_i =
	      Aztec.update_index[Ncomp + Nrho_bar +
				 Nunk_per_node * loc_inode];
	  /* printf("inode: %d icomp: %d  elec_pot: %9.6f \n",inode,Ncomp,x[loc_i]); */
/*          printf(" %9.6f  %9.6f  %9.6f \n",nodepos[0],nodepos[1],x[loc_i]);*/
/*          printf(" %9.6f  %9.6f  ",nodepos[0],nodepos[1]);*/

	  /*for (iel_w=0; iel_w<Nelems_S[Nlists_HW-1][loc_inode]; iel_w++){
	     surf_norm = Surf_normal[Nlists_HW-1][loc_inode][iel_w]; */
	  for (iel_w = 0; iel_w < Nelems_S[0][loc_inode]; iel_w++)
	    {
	      surf_norm = Surf_normal[0][loc_inode][iel_w];
	      el_type = Surf_elem_type[loc_inode][iel_w];
/*printf("el_type %d  %d %d\n", el_type ,loc_inode,iel_w);*/
	      idim = abs (surf_norm) - 1;

	      for (jdim = 0; jdim < Ndim; jdim++)
		{
		  if (jdim == idim)
		    {
		      if (surf_norm < 0)
			deriv_x[iel_w][jdim] =
			  calc_deriv (idim, inode_box, BFD, &blocked, x, 0);
		      else
			deriv_x[iel_w][jdim] =
			  calc_deriv (idim, inode_box, FFD, &blocked, x, 0);
/*                    printf("%d  %9.6f",jdim,deriv_x[iel_w][jdim]);*/
		    }
		  else
		    {
		      find_offset (el_type, jdim, offset);
		      jnode_box =
			offset_to_node_box (ijk_box, offset, reflect_flag);
		      /*if (Nodes_2_boundary_wall[Nlists_HW-1][jnode_box] == -1){ */
		      if (Nodes_2_boundary_wall[0][jnode_box] == -1)
			{
			  /*     printf("trouble ... the derivatives are not within surface elements !!"); */
			}
		      if (Matrix_fill_flag < 3)
			loc_j =
			  B2L_unknowns[jnode_box * Nunk_per_node + Ncomp];
		      else
			loc_j =
			  B2L_unknowns[jnode_box * Nunk_per_node + Ncomp +
				       Nrho_bar];

		      deriv_x[iel_w][jdim] =
			offset[jdim] * (x[loc_j] - x[loc_i]) / Esize_x[jdim];
/*                    printf("%d  %9.6f",jdim,deriv_x[iel_w][jdim]);*/
		    }

		  dot_prod[iel_w] +=
		    deriv_x[iel_w][jdim] * deriv_x[iel_w][jdim];
		}
	    }
/*          printf("\n");*/


	  /*for (iel_w=0; iel_w<Nelems_S[Nlists_HW-1][loc_inode]; iel_w++){
	     surf_norm = Surf_normal[Nlists_HW-1][loc_inode][iel_w]; */

	  for (iel_w = 0; iel_w < Nelems_S[0][loc_inode]; iel_w++)
	    {

	      jwall = Surf_elem_to_wall[0][loc_inode][iel_w];
	      surf_norm = Surf_normal[0][loc_inode][iel_w];
	      idim = abs (surf_norm) - 1;
	      prefac = (double) (surf_norm / abs (surf_norm)) *
		Area_surf_el[idim] / Nnodes_per_el_S;

	      integrand =
		deriv_x[iel_w][idim] * deriv_x[iel_w][idim] -
		0.5 * dot_prod[iel_w];

	      Sum_dphi_dx[jwall][idim] -= prefac * integrand;
	      if (idim == 1)
		store1 += prefac * integrand;

	      for (jdim = 0; jdim < Ndim; jdim++)
		{

		  if (idim != jdim)
		    integrand = deriv_x[iel_w][jdim] * deriv_x[iel_w][idim];
		  else
		    integrand = 0.0;

		  Sum_dphi_dx[jwall][jdim] -= prefac * integrand;
		  if (jdim == 1)
		    store2 -= prefac * integrand;
		}


	    }			/* end of surface element loop */
	}			/* end of test for boundary node */
    }				/* end of loop over nodes on this processor */

  store1_tot = AZ_gsum_double (store1,
			       Aztec.proc_config) * Temp_elec / (4.0 * PI *
								 0.05);
  store2_tot =
    AZ_gsum_double (store2,
		    Aztec.proc_config) * Temp_elec / (4.0 * PI * 0.05);
  if (Proc == 0)
    printf ("normal contribution: %9.6f  tangential contribution:%9.6f",
	    store1_tot, store2_tot);

  return;
}

/*************************************************************************
find_offset: in this routine we use the current node, the surface element
             normal, and the surface element type to determine the offset
             for the appropriate derivative */

void
find_offset (int el_type, int jdim, int *offset)
{
  offset[0] = offset[1] = offset[2] = 0;

  if (jdim == 0)
    {
      if (el_type == RIGHT_BACK ||
	  el_type == RIGHT_FRONT ||
	  el_type == RIGHT_UP || el_type == RIGHT_DOWN)
	offset[0] = 1;

      else if (el_type == LEFT_BACK ||
	       el_type == LEFT_FRONT ||
	       el_type == LEFT_UP || el_type == LEFT_DOWN)
	offset[0] = -1;
    }
  else if (jdim == 1)
    {
      if (el_type == UP_BACK ||
	  el_type == UP_FRONT || el_type == RIGHT_UP || el_type == LEFT_UP)
	offset[1] = 1;

      else if (el_type == DOWN_BACK ||
	       el_type == DOWN_FRONT ||
	       el_type == RIGHT_DOWN || el_type == LEFT_DOWN)
	offset[1] = -1;
    }
  else if (jdim == 2)
    {
      if (el_type == UP_FRONT ||
	  el_type == RIGHT_FRONT ||
	  el_type == LEFT_FRONT || el_type == DOWN_FRONT)
	offset[2] = 1;

      if (el_type == UP_BACK ||
	  el_type == RIGHT_BACK ||
	  el_type == LEFT_BACK || el_type == DOWN_BACK)
	offset[2] = -1;
    }
  return;
}

/*********************************************************************
 * calc_deriv : calculate a derivative of the electric potential!!   */

double
calc_deriv (int idim, int inode0, int flag, int *blocked, double *x,
	    int ilist)
{
  int inode1, inode2, offset1[3], offset2[3],
    loc_i0, loc_i1, loc_i2, iwall1, iwall2, jdim, ijk_box[3], reflect_flag[3];
  double deriv = 0.0;

  node_box_to_ijk_box (inode0, ijk_box);

  for (jdim = 0; jdim < Ndim; jdim++)
    {
      offset1[jdim] = 0;
      offset2[jdim] = 0;
    }

  switch (flag)
    {
    case CFD:
      offset1[idim] = -1;
      offset2[idim] = 1;
      break;

    case FFD:
      offset1[idim] = 1;
      offset2[idim] = 2;
      break;

    case BFD:
      offset1[idim] = -1;
      offset2[idim] = -2;
      break;
    }

  inode1 = offset_to_node_box (ijk_box, offset1, reflect_flag);
  inode2 = offset_to_node_box (ijk_box, offset2, reflect_flag);

  if (Matrix_fill_flag < 3)
    {
      loc_i0 = B2L_unknowns[inode0 * Nunk_per_node + Ncomp];
      loc_i1 = B2L_unknowns[inode1 * Nunk_per_node + Ncomp];
      loc_i2 = B2L_unknowns[inode2 * Nunk_per_node + Ncomp];
    }
  else
    {
      loc_i0 = B2L_unknowns[inode0 * Nunk_per_node + Ncomp + Nrho_bar];
      loc_i1 = B2L_unknowns[inode1 * Nunk_per_node + Ncomp + Nrho_bar];
      loc_i2 = B2L_unknowns[inode2 * Nunk_per_node + Ncomp + Nrho_bar];
    }

/*   iwall1 = Nodes_2_boundary_wall[Nlists_HW-1][inode1];
   iwall2 = Nodes_2_boundary_wall[Nlists_HW-1][inode2];*/

  iwall1 = Nodes_2_boundary_wall[ilist][inode1];
  iwall2 = Nodes_2_boundary_wall[ilist][inode2];

  if (iwall1 == -1 && iwall2 == -1)
    {
      *blocked = FALSE;
      switch (flag)
	{
	case CFD:
	  deriv = x[loc_i2] - x[loc_i1];
	  break;
	case FFD:
	  deriv = -3 * x[loc_i0] + 4 * x[loc_i1] - x[loc_i2];
	  break;
	case BFD:
	  deriv = 3 * x[loc_i0] - 4 * x[loc_i1] + x[loc_i2];
	  break;
	}
      deriv /= (2.0 * Esize_x[idim]);
    }
  else
    {
      *blocked = TRUE;
      switch (flag)
	{
	case FFD:
	  deriv = x[loc_i1] - x[loc_i0];
	  break;
	case BFD:
	  deriv = x[loc_i0] - x[loc_i1];
	  break;
	}
      deriv /= (Esize_x[idim]);
    }

  return deriv;
}

/***************************************************************************
 * find_pot_derivs: For the hard wall case, sum the densities at the wall to  *
 *              get the force on the wall.                                 */

void
find_pot_derivs (double *x, double *psi_deriv)
{
  int loc_i, loc_inode, iwall, idim, iunk, ilist, iel_w, inode, surf_norm, i;
  double prefac;

  if (Matrix_fill_flag < 3)
    iunk = Ncomp;
  else
    iunk = Ncomp + Nrho_bar;

  ilist = 0;
  psi_deriv[0] = 0.0;
  psi_deriv[1] = 0.0;

  for (loc_inode = 0; loc_inode < Nnodes_per_proc; loc_inode++)
    {
      inode = Aztec.update[Nunk_per_node * loc_inode] / Nunk_per_node;
      iwall = Nodes_2_boundary_wall[ilist][node_to_node_box (inode)];

      if (iwall != -1)
	{

	  loc_i = Aztec.update_index[iunk + Nunk_per_node * loc_inode];
	  for (iel_w = 0; iel_w < Nelems_S[ilist][loc_inode]; iel_w++)
	    {

	      surf_norm = Surf_normal[ilist][loc_inode][iel_w];
	      idim = abs (surf_norm) - 1;

	      prefac =
		(double) (surf_norm / abs (surf_norm)) * Area_surf_el[idim];

	      if (prefac > 0.0)
		{		/* RHS */
		  i = 1;
		  psi_deriv[i] = (-x[iunk + Nunk_per_node * (loc_inode + 2)]
				  + 4. * x[iunk +
					   Nunk_per_node * (loc_inode + 1)] -
				  3. * x[loc_i]) / (2. * Esize_x[0]);
		}
	      else
		{		/* LHS */
		  i = 0;
		  psi_deriv[i] = (3. * x[loc_i]
				  - 4. * x[iunk +
					   Nunk_per_node * (loc_inode - 1)] +
				  x[iunk +
				    Nunk_per_node * (loc_inode -
						     2)]) / (2.0 *
							     Esize_x[0]);
		}


	    }			/* end of surface element loop */
	}			/* end of test for boundary node */
    }				/* end of loop over nodes on this processor */
  return;
}

/***************************************************************************
 * sum_rho_midplane: For the hard wall case, sum the densities at the wall to  *
 *                   get the force on the wall.                                */

double
sum_rho_midplane (double *x)
{
  int loc_i, loc_inode, icomp, inode, ijk[3];
  double rho_mid_sum;

  rho_mid_sum = 0.0;

  for (loc_inode = 0; loc_inode < Nnodes_per_proc; loc_inode++)
    {
      inode = Aztec.update[Nunk_per_node * loc_inode] / Nunk_per_node;
      node_to_ijk (inode, ijk);

      if (ijk[0] == Nodes_x[0] - 1)
	{
	  icomp = 0;
	  loc_i = Aztec.update_index[icomp + Nunk_per_node * loc_inode];
	  if (Ipot_ff_n == IDEAL_GAS)
	    {
	      for (icomp = 0; icomp < Ncomp; icomp++)
		{
		  loc_i =
		    Aztec.update_index[icomp + Nunk_per_node * loc_inode];
		  if (Ndim == 1)
		    rho_mid_sum += x[loc_i];
		  else if (Ndim == 2)
		    rho_mid_sum += x[loc_i] * Esize_x[1];
		}
	    }
	  else
	    {
	      rho_mid_sum += calc_local_pressure (&(x[loc_i]));
	    }

	}			/* end of loop over nodes on this processor */
    }				/* end of icomp loop */
  return (rho_mid_sum);
}

/****************************************************************************
calc_geom_factor:  In this routine calculate the geometrical factor associated
                   with the surfaces of interest.
void calc_geom_factor(int iwall, double *nodepos,double *factor)
{
  double r_iw[3],r_iw_mag,r_iw_mag_sq;
  int idim;

   * assume that this node is point 3, the center of surface # 1 is 1,
   * the center of surface #2 is 2.  Calculate the magnitudes of the vectors
   * R_12 and r_23. 
   * 

   if (Ndim == 1 || Surface_type[0] == 0) 
      for (idim=0; idim<Ndim; idim++) factor[0] = 1.0;

   else {
     for (idim=0; idim<Ndim; idim++) 
          r_iw[idim] = WallPos[idim][iwall]-nodepos[idim];
  
     r_iw_mag_sq = 0.0;
     for (idim=0; idim<Ndim; idim++) 
           r_iw_mag_sq += r_iw[idim]*r_iw[idim];

     r_iw_mag = sqrt(r_iw_mag_sq);

     for (idim=0; idim<Ndim; idim++) 
           factor[idim] = fabs(r_iw[idim])/r_iw_mag;
   }
  return;
}
*/
/****************************************************************************
*integrate_rho_vdash: Find p_tilde as the integral of rho_vdash 
*                        for the 1-dimensional problem.                  */
void
integrate_rho_vdash (double *x, double **rho_vdash)
{
  int loc_i, loc_inode, iwall, idim, icomp, ilist,
    iunk, inode, reflect_flag[3], ijk[3], iel, ielement, nel_hit;

  for (iwall = 0; iwall < Nwall; iwall++)
    {
      for (idim = 0; idim < Ndim; idim++)
	{
	  rho_vdash[iwall][idim] = 0.0;
	}
    }

  for (idim = 0; idim < Ndim; idim++)
    reflect_flag[idim] = FALSE;

  for (loc_inode = 0; loc_inode < Nnodes_per_proc; loc_inode++)
    {

      inode = Aztec.update[Nunk_per_node * loc_inode] / Nunk_per_node;
      node_to_ijk (inode, ijk);

      for (icomp = 0; icomp < Ncomp; icomp++)
	{
	  for (iwall = 0; iwall < Nwall; iwall++)
	    {


	      ilist = 0;
	      nel_hit = Nnodes_per_el_V;
	      for (iel = 0; iel < Nnodes_per_el_V; iel++)
		{
		  ielement = node_to_elem (inode, iel, reflect_flag);
		  if (ielement != -2)
		    {
		      if (ielement == -1)
			nel_hit--;
		      else if (Wall_elems[ilist][ielement] != -1)
			nel_hit--;
		    }
		}

	      for (idim = 0; idim < Ndim; idim++)
		{
		  if (ijk[idim] == 0 && (Type_bc[idim][0] == REFLECT ||
					 Type_bc[idim][0] == IN_BULK ||
					 Type_bc[idim][0] == LAST_NODE))
		    nel_hit /= 2;
		  else if (ijk[idim] == Nodes_x[idim] - 1 &&
			   (Type_bc[idim][1] == REFLECT ||
			    Type_bc[idim][1] == IN_BULK ||
			    Type_bc[idim][1] == LAST_NODE))
		    nel_hit /= 2;

		  loc_i =
		    Aztec.update_index[icomp + Nunk_per_node * loc_inode];
		  iunk = iwall * Ncomp + icomp;
		  rho_vdash[iwall][idim] -=
		    (x[loc_i] * Vext_dash[loc_inode][iunk][idim])
		    * nel_hit * Vol_el / ((double) Nnodes_per_el_V);
		}		/* end of idim loop */
	    }			/* end of Nwall loop */

	}			/* end of icomp loop */
    }				/* end of loc_inode loop */
}

/****************************************************************************/
/*calc_local_pressure:  This routine calculates the local pressure in the 
                      solution */

double
calc_local_pressure (double *rho)
{
  int icomp;
  double pi6, hs_diam_cubed, xsi0, xsi1, xsi2, xsi3, y1, y2, y3, betap_hs;

  xsi0 = 0.0;
  xsi1 = 0.0;
  xsi2 = 0.0;
  xsi3 = 0.0;
  pi6 = PI / 6.0;		/* shorthand  for pi/6                 */

  /*  Determine the effective hard sphere diameters 
     For now we will set these to unity, but in the future we
     can define a temperature dependent diameter. Doing so
     offers a way to provide a better mean field equation of
     state. In essence, for  an attractive (e.g LJ) fluid (mixture)
     we can use the effective hard sphere diameter to off set the 
     shortcomings of the PY + mean field approximation. See the
     work by Telo da Gama et al.. 
   */

  /* calculate the hard sphere diamtere ... this can be 
     turned into a T-dependent diameter */

  for (icomp = 0; icomp < Ncomp; ++icomp)
    Hs_diam[icomp] = 1.0;

  /*  calculate the constants xsi and introduce some shorthand */

  for (icomp = 0; icomp < Ncomp; ++icomp)
    {
      hs_diam_cubed = POW_DOUBLE_INT (Hs_diam[icomp], 3);
      xsi0 += pi6 * rho[icomp] * hs_diam_cubed;
      xsi1 += pi6 * rho[icomp] * hs_diam_cubed * Sigma_ff[icomp][icomp];
      xsi2 +=
	pi6 * rho[icomp] * hs_diam_cubed *
	POW_DOUBLE_INT (Sigma_ff[icomp][icomp], 2);
      xsi3 +=
	pi6 * rho[icomp] * hs_diam_cubed *
	POW_DOUBLE_INT (Sigma_ff[icomp][icomp], 3);
    }
  y1 = 1.0 - xsi3;
  y2 = y1 * y1;
  y3 = y1 * y1 * y1;

  /* the hard sphere pressure in units of kT and Sigma_ff[1]^3 */

  betap_hs = (1.0 / pi6) * (xsi0 / y1 + 3.0 * xsi1 * xsi2 / y2 +
			    3.0 * POW_DOUBLE_INT (xsi2, 3) / y3);
  return (betap_hs);
}

/***********************************************************************/
