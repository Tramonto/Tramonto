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
 *  FILE: dft_fill_rosen.c
 *
 *  This file contains the fill routines for the rho-only
 *  rosenfeld formulation of the equations.
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"

void load_Jacobian (double *, struct RB_Struct, int, int, int *, int, int);
void load_Jac_Save_Jacobian (double *, struct RB_Struct, int, int,
			     int *, int, double ***, int ***);

void fast_fill_Jacobian (int, int,
			 int, int, double, double *, double *, int *, int);
void medium_fill_Jacobian (int, int, int,
			   int, int, double, double *, double *, int *, int);
void slow_fill_Jacobian (int, int, int, int, int,
			 double, double *, double *, int *, int);

void final_Jacobian_fill (int, int, int, double, int *, double,
			  double *, double *, int *, int);

struct RB_Struct d2phi_drb2_delta (struct RB_Struct, double, int *, double,
				   double, double, double *, int);
struct RB_Struct d2phi_drb2_delta2 (struct RB_Struct, double, int *, double,
				    double, double, double *, int);
struct RB_Struct d2phi_drb2_theta (struct RB_Struct, double, int *, int);
struct RB_Struct d2phi_drb2_theta2 (struct RB_Struct, double, int *, int);

/****************************************************************************/
/* load_nonlocal_hs_rosen: Here we load all the dphi_drb terms for the 
                        Rosenfeld functional.                         */
double
load_nonlocal_hs_rosen (int sten_type, int loc_i, int icomp,
			int izone, int *ijk_box, int fill_flag,
			struct RB_Struct *rho_bar,
			struct RB_Struct *dphi_drb,
			struct RB_Struct *dphi_drb_bulk,
			struct RB_Struct *dphi_drb_bulk_left,
			struct RB_Struct *dphi_drb_bulk_right,
			double *resid, double *mat_row, int *bindx_tmp,
			double ***jac_weights_hs, int ***jac_columns_hs,
			int resid_only_flag)
{
  int **sten_offset, *offset, isten;
  double *sten_weight, weight;
  struct Stencil_Struct *sten;

  int jzone = 0, loc_jnode, jnode_box, idim, jzone_mesh;
  int reflect_flag[NDIM_MAX];
  double sign[3], t_lj = 0.0;
  struct RB_Struct tmp;
  double resid_start;
  if (fill_flag != MSR_PREPROCESS)
    resid_start = resid[loc_i];

  tmp.S0 = tmp.S1 = tmp.S2 = tmp.S3 = tmp.V1[0] = tmp.V1[1]
    = tmp.V1[2] = tmp.V2[0] = tmp.V2[1] = tmp.V2[2] = 0.0;

  sten = &(Stencil[sten_type][izone][icomp]);
  sten_offset = sten->Offset;
  sten_weight = sten->Weight;

  for (isten = 0; isten < sten->Length; isten++)
    {

      /* get offset and weight for this stencil point */
      offset = sten_offset[isten];
      weight = sten_weight[isten];

      /* Find the Stencil point */
      jnode_box = offset_to_node_box (ijk_box, offset, reflect_flag);

      if (jnode_box >= 0)
	{

	  /* set zone based on position in domain and Nzone from input file */

	  jzone_mesh = Nodes_to_zone[jnode_box];
	  if (Nzone > 1 && Coarser_jac > 0 && jzone_mesh < Nzone - 1)
	    {
	      if (Coarser_jac == 1)
		{
		  if (jzone_mesh == 0)
		    jzone = jzone_mesh + 1;
		  else
		    jzone = jzone_mesh;
		}
	      else if (Coarser_jac == 2)
		jzone = jzone_mesh + 1;
	      else if (Coarser_jac == 3)
		jzone = Nzone - 1;
	      else if (Coarser_jac == 4)
		jzone = Nzone - 2;
	      else if (Coarser_jac == 5)
		jzone = Nzone - 1;
	    }
	  else
	    jzone = jzone_mesh;

	  if (fill_flag != MSR_PREPROCESS)
	    {
	      loc_jnode = B2L_1stencil[jnode_box];

	      if (sten_type == DELTA_FN)
		{
		  resid[loc_i] += weight *
		    (dphi_drb[loc_jnode].S0 * Inv_4pirsq[icomp] +
		     dphi_drb[loc_jnode].S1 * Inv_4pir[icomp] +
		     dphi_drb[loc_jnode].S2);
		  for (idim = 0; idim < Ndim; idim++)
		    {
		      if (reflect_flag[idim] == FALSE)
			sign[idim] = 1.0;
		      else
			sign[idim] = -1.0;
		      resid[loc_i] -= sign[idim] * weight *
			(dphi_drb[loc_jnode].V1[idim] * Inv_4pir[icomp]
			 + dphi_drb[loc_jnode].V2[idim]) *
			(offset[idim] * Esize_x[idim] * Inv_rad[icomp]);
		    }
		  if (!resid_only_flag)
		    {
		      if (Type_func == 0)
			tmp =
			  d2phi_drb2_delta (rho_bar[loc_jnode], weight,
					    offset, Inv_rad[icomp],
					    Inv_4pir[icomp],
					    Inv_4pirsq[icomp], sign,
					    fill_flag);
		      else
			tmp =
			  d2phi_drb2_delta2 (rho_bar[loc_jnode], weight,
					     offset, Inv_rad[icomp],
					     Inv_4pir[icomp],
					     Inv_4pirsq[icomp], sign,
					     fill_flag);
		    }
		}
	      else if (sten_type == THETA_FN)
		{
		  resid[loc_i] += weight * dphi_drb[loc_jnode].S3;
		  if (!resid_only_flag)
		    {
		      if (Type_func == 0)
			tmp =
			  d2phi_drb2_theta (rho_bar[loc_jnode], weight,
					    offset, fill_flag);
		      else
			tmp =
			  d2phi_drb2_theta2 (rho_bar[loc_jnode], weight,
					     offset, fill_flag);
		    }
		}

	    }			/* end of  if (fill_flag != MSR_PREPROCESS) loop */

	  t_lj -= MPI_Wtime ();
	  if (!resid_only_flag)
	    {
	      if (Matrix_fill_flag == FULL_MSR_FILL)
		load_Jacobian (mat_row, tmp, jnode_box, jzone, bindx_tmp,
			       Fast_fill_TF[jnode_box], fill_flag);
	      else
		load_Jac_Save_Jacobian (mat_row, tmp, jnode_box, jzone,
					bindx_tmp, fill_flag, jac_weights_hs,
					jac_columns_hs);
	    }
	  t_lj += MPI_Wtime ();
	}
      else if (fill_flag != MSR_PREPROCESS && (jnode_box == -1 ||
					       jnode_box == -3 ||
					       jnode_box == -4))
	{
	  if (jnode_box == -1)
	    {
	      if (sten_type == DELTA_FN)
		resid[loc_i] += weight *
		  (dphi_drb_bulk->S0 * Inv_4pirsq[icomp] +
		   dphi_drb_bulk->S1 * Inv_4pir[icomp] + dphi_drb_bulk->S2);
	      else if (sten_type == THETA_FN)
		{
		  resid[loc_i] += weight * dphi_drb_bulk->S3;
		}
	    }
	  else if (jnode_box == -3)
	    {
	      if (sten_type == DELTA_FN)
		resid[loc_i] += weight *
		  (dphi_drb_bulk_left->S0 * Inv_4pirsq[icomp] +
		   dphi_drb_bulk_left->S1 * Inv_4pir[icomp] +
		   dphi_drb_bulk_left->S2);
	      else if (sten_type == THETA_FN)
		resid[loc_i] += weight * dphi_drb_bulk_left->S3;
	    }
	  if (jnode_box == -4)
	    {
	      if (sten_type == DELTA_FN)
		resid[loc_i] += weight *
		  (dphi_drb_bulk_right->S0 * Inv_4pirsq[icomp] +
		   dphi_drb_bulk_right->S1 * Inv_4pir[icomp] +
		   dphi_drb_bulk_right->S2);
	      else if (sten_type == THETA_FN)
		resid[loc_i] += weight * dphi_drb_bulk_right->S3;
	    }
	}
    }
  return t_lj;
}

/*****************************************************************************/
/* pre_calc_dphi_drb: rho_bars are calculated at each wall & fluid node*/

void
pre_calc_dphi_drb (struct RB_Struct *dphi_drb,
		   struct RB_Struct *rho_bar,
		   struct RB_Struct *dphi_drb_bulk_ptr,
		   struct RB_Struct *dphi_drb_bulk_left_ptr,
		   struct RB_Struct *dphi_drb_bulk_right_ptr)
{
  double rb0, rb1, rb2, rb3, rb1v[3], rb2v[3], inv_one_m_rb3,
    inv_one_m_rb3_sq, inv_one_m_rb3_3rd, DOT_rho12, DOT_rho22;
  double rb0_l = 0.0, rb1_l = 0.0, rb2_l = 0.0, rb3_l = 0.0;
  double rb0_r = 0.0, rb1_r = 0.0, rb2_r = 0.0, rb3_r = 0.0;
  int idim, icomp, loc_inode;
  int inode_box;
  int i, imax;
  double rho;

  for (inode_box = 0; inode_box < Nnodes_box; inode_box++)
    {
      if (B2L_1stencil[inode_box] != -1)
	{
	  loc_inode = B2L_1stencil[inode_box];
	  rb3 = rho_bar[loc_inode].S3;
	  rb2 = rho_bar[loc_inode].S2;
	  rb1 = rho_bar[loc_inode].S1;
	  rb0 = rho_bar[loc_inode].S0;

	  for (idim = 0; idim < Ndim; idim++)
	    {
	      rb2v[idim] = rho_bar[loc_inode].V2[idim];
	      rb1v[idim] = rho_bar[loc_inode].V1[idim];
	    }
	  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
	  inv_one_m_rb3_sq = inv_one_m_rb3 * inv_one_m_rb3;
	  inv_one_m_rb3_3rd = inv_one_m_rb3_sq * inv_one_m_rb3;

	  dphi_drb[loc_inode].S0 = log (inv_one_m_rb3);
	  dphi_drb[loc_inode].S1 = rb2 * inv_one_m_rb3;
	  dphi_drb[loc_inode].S2 = rb1 * inv_one_m_rb3 +
	    rb2 * rb2 * inv_one_m_rb3_sq / (8.0 * PI);
	  dphi_drb[loc_inode].S3 = rb0 * inv_one_m_rb3 +
	    rb1 * rb2 * inv_one_m_rb3_sq +
	    rb2 * rb2 * rb2 * inv_one_m_rb3_3rd / (12.0 * PI);
	  DOT_rho12 = 0.0;
	  DOT_rho22 = 0.0;
	  for (idim = 0; idim < Ndim; idim++)
	    {
	      DOT_rho12 += rb1v[idim] * rb2v[idim];
	      DOT_rho22 += rb2v[idim] * rb2v[idim];

	      dphi_drb[loc_inode].V1[idim] = -rb2v[idim] * inv_one_m_rb3;
	      dphi_drb[loc_inode].V2[idim] = -rb1v[idim] * inv_one_m_rb3 -
		rb2 * rb2v[idim] * inv_one_m_rb3_sq / (4.0 * PI);
	    }

	  dphi_drb[loc_inode].S2 +=
	    -DOT_rho22 * inv_one_m_rb3_sq / (8.0 * PI);
	  dphi_drb[loc_inode].S3 +=
	    -DOT_rho12 * inv_one_m_rb3_sq -
	    rb2 * DOT_rho22 * inv_one_m_rb3_3rd / (4.0 * PI);
	}			/* End B2L if statement */
    }

  /* Calculate dphi_drb_bulk, values in the bulk */

  rb0 = rb1 = rb2 = rb3 = 0.0;
  for (icomp = 0; icomp < Ncomp; icomp++)
    {
      rb0 += Rho_b[icomp];
      rb1 += Rho_b[icomp] * Sigma_ff[icomp][icomp] / 2.0;
      rb2 +=
	Rho_b[icomp] * PI * Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
      rb3 +=
	Rho_b[icomp] * (PI / 6.0) * Sigma_ff[icomp][icomp] *
	Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
    }
  if (Iliq_vap == 3 || Lsteady_state == TRUE)
    {
      if (Lsteady_state == TRUE)
	imax = Ncomp;
      else
	imax = 1;

      rb0_l = rb1_l = rb2_l = rb3_l = 0.0;
      rb0_r = rb1_r = rb2_r = rb3_r = 0.0;

      for (i = 0; i < imax; i++)
	{
	  if (Lsteady_state == TRUE)
	    rho = Rho_b_LBB[i];
	  else
	    rho = Rho_coex[1];

	  rb0_l += rho;
	  rb1_l += rho * Sigma_ff[i][i] / 2.0;
	  rb2_l += rho * PI * Sigma_ff[i][i] * Sigma_ff[i][i];
	  rb3_l += rho * (PI / 6.0) * Sigma_ff[i][i]
	    * Sigma_ff[i][i] * Sigma_ff[i][i];

	  if (Lsteady_state == TRUE)
	    rho = Rho_b_RTF[i];
	  else
	    rho = Rho_coex[0];

	  rb0_r += rho;
	  rb1_r += rho * Sigma_ff[i][i] / 2.0;
	  rb2_r += rho * PI * Sigma_ff[i][i] * Sigma_ff[i][i];
	  rb3_r += rho * (PI / 6.0) * Sigma_ff[i][i]
	    * Sigma_ff[i][i] * Sigma_ff[i][i];
	}
    }

  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3 * inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq * inv_one_m_rb3;

  dphi_drb_bulk_ptr->S0 = log (inv_one_m_rb3);
  dphi_drb_bulk_ptr->S1 = rb2 * inv_one_m_rb3;
  dphi_drb_bulk_ptr->S2 = rb1 * inv_one_m_rb3 +
    rb2 * rb2 * inv_one_m_rb3_sq / (8.0 * PI);
  dphi_drb_bulk_ptr->S3 = rb0 * inv_one_m_rb3 +
    rb1 * rb2 * inv_one_m_rb3_sq +
    rb2 * rb2 * rb2 * inv_one_m_rb3_3rd / (12.0 * PI);

  /* dphi_drb_bulk_ptr->V1 = dphi_drb_bulk_ptr->V2 = 0 in the bulk */
  /* -- not calculated */

  if (Iliq_vap == 3 || Lsteady_state)
    {
      inv_one_m_rb3 = 1.0 / (1.0 - rb3_l);
      inv_one_m_rb3_sq = inv_one_m_rb3 * inv_one_m_rb3;
      inv_one_m_rb3_3rd = inv_one_m_rb3_sq * inv_one_m_rb3;

      dphi_drb_bulk_left_ptr->S0 = log (inv_one_m_rb3);
      dphi_drb_bulk_left_ptr->S1 = rb2_l * inv_one_m_rb3;
      dphi_drb_bulk_left_ptr->S2 = rb1_l * inv_one_m_rb3 +
	rb2_l * rb2_l * inv_one_m_rb3_sq / (8.0 * PI);
      dphi_drb_bulk_left_ptr->S3 = rb0_l * inv_one_m_rb3 +
	rb1_l * rb2_l * inv_one_m_rb3_sq +
	rb2_l * rb2_l * rb2_l * inv_one_m_rb3_3rd / (12.0 * PI);

      inv_one_m_rb3 = 1.0 / (1.0 - rb3_r);
      inv_one_m_rb3_sq = inv_one_m_rb3 * inv_one_m_rb3;
      inv_one_m_rb3_3rd = inv_one_m_rb3_sq * inv_one_m_rb3;

      dphi_drb_bulk_right_ptr->S0 = log (inv_one_m_rb3);
      dphi_drb_bulk_right_ptr->S1 = rb2_r * inv_one_m_rb3;
      dphi_drb_bulk_right_ptr->S2 = rb1_r * inv_one_m_rb3 +
	rb2_r * rb2_r * inv_one_m_rb3_sq / (8.0 * PI);
      dphi_drb_bulk_right_ptr->S3 = rb0_r * inv_one_m_rb3 +
	rb1_r * rb2_r * inv_one_m_rb3_sq +
	rb2_r * rb2_r * rb2_r * inv_one_m_rb3_3rd / (12.0 * PI);
    }

}

/*****************************************************************************/
/* pre_calc_dphi_drb2: rho_bars are calculated at each wall & fluid node 
                       using new Rosenfeld functionals !                     */

void
pre_calc_dphi_drb2 (struct RB_Struct *dphi_drb,
		    struct RB_Struct *rho_bar,
		    struct RB_Struct *dphi_drb_bulk_ptr,
		    struct RB_Struct *dphi_drb_bulk_left_ptr,
		    struct RB_Struct *dphi_drb_bulk_right_ptr)
{
  double rb0, rb1, rb2, rb3, rb1v[3], rb2v[3], inv_one_m_rb3,
    inv_one_m_rb3_sq, inv_one_m_rb3_3rd, DOT_rho12, DOT_rho22;
  double rb0_l = 0.0, rb1_l = 0.0, rb2_l = 0.0, rb3_l = 0.0;
  double rb0_r = 0.0, rb1_r = 0.0, rb2_r = 0.0, rb3_r = 0.0;
  int idim, icomp, loc_inode;
  int inode_box;
  int i, imax;
  double rho, alpha, alpha_sq, alpha_cb, beta, gamma[3];

  for (inode_box = 0; inode_box < Nnodes_box; inode_box++)
    {
      if (B2L_1stencil[inode_box] != -1)
	{
	  loc_inode = B2L_1stencil[inode_box];
	  rb3 = rho_bar[loc_inode].S3;
	  rb2 = rho_bar[loc_inode].S2;
	  rb1 = rho_bar[loc_inode].S1;
	  rb0 = rho_bar[loc_inode].S0;

	  for (idim = 0; idim < Ndim; idim++)
	    {
	      rb2v[idim] = rho_bar[loc_inode].V2[idim];
	      rb1v[idim] = rho_bar[loc_inode].V1[idim];
	    }
	  DOT_rho12 = 0.0;
	  DOT_rho22 = 0.0;
	  for (idim = 0; idim < Ndim; idim++)
	    {
	      DOT_rho12 += rb1v[idim] * rb2v[idim];
	      DOT_rho22 += rb2v[idim] * rb2v[idim];
	    }
	  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
	  inv_one_m_rb3_sq = inv_one_m_rb3 * inv_one_m_rb3;
	  inv_one_m_rb3_3rd = inv_one_m_rb3_sq * inv_one_m_rb3;

	  /* same as all old rosenfeld functional contributions */
	  dphi_drb[loc_inode].S0 = log (inv_one_m_rb3);
	  dphi_drb[loc_inode].S1 = rb2 * inv_one_m_rb3;
	  dphi_drb[loc_inode].S2 = rb1 * inv_one_m_rb3;
	  dphi_drb[loc_inode].S3 = rb0 * inv_one_m_rb3 +
	    (rb1 * rb2 - DOT_rho12) * inv_one_m_rb3_sq;

	  for (idim = 0; idim < Ndim; idim++)
	    {
	      dphi_drb[loc_inode].V1[idim] = -rb2v[idim] * inv_one_m_rb3;
	      dphi_drb[loc_inode].V2[idim] = -rb1v[idim] * inv_one_m_rb3;
	    }

	  /* new rosenfeld functional contributions */
	  if (rb2 > 1.e-15)
	    {
	      alpha = rb2 - DOT_rho22 / rb2;
	      beta = 1.0 + DOT_rho22 / (rb2 * rb2);
	      for (idim = 0; idim < Ndim; idim++)
		{
		  gamma[idim] = rb2v[idim] / rb2;
		}
	    }
	  else
	    {
	      alpha = rb2;
	      beta = 1.0;
	      for (idim = 0; idim < Ndim; idim++)
		gamma[idim] = 0.0;
	    }
	  alpha_sq = alpha * alpha;
	  alpha_cb = alpha_sq * alpha;

	  dphi_drb[loc_inode].S2 +=
	    alpha_sq * beta * inv_one_m_rb3_sq / (8.0 * PI);
	  dphi_drb[loc_inode].S3 +=
	    alpha_cb * inv_one_m_rb3_3rd / (12.0 * PI);

	  for (idim = 0; idim < Ndim; idim++)
	    dphi_drb[loc_inode].V2[idim] -=
	      inv_one_m_rb3_sq * alpha_sq * gamma[idim] / (4.0 * PI);

	}			/* End B2L if statement */
    }


  /* Calculate dphi_drb_bulk, values in the bulk */

  rb0 = rb1 = rb2 = rb3 = 0.0;
  for (icomp = 0; icomp < Ncomp; icomp++)
    {
      rb0 += Rho_b[icomp];
      rb1 += Rho_b[icomp] * Sigma_ff[icomp][icomp] / 2.0;
      rb2 +=
	Rho_b[icomp] * PI * Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
      rb3 +=
	Rho_b[icomp] * (PI / 6.0) * Sigma_ff[icomp][icomp] *
	Sigma_ff[icomp][icomp] * Sigma_ff[icomp][icomp];
    }
  if (Iliq_vap == 3 || Lsteady_state == TRUE)
    {
      if (Lsteady_state == TRUE)
	imax = Ncomp;
      else
	imax = 1;

      rb0_l = rb1_l = rb2_l = rb3_l = 0.0;
      rb0_r = rb1_r = rb2_r = rb3_r = 0.0;

      for (i = 0; i < imax; i++)
	{
	  if (Lsteady_state == TRUE)
	    rho = Rho_b_LBB[i];
	  else
	    rho = Rho_coex[1];

	  rb0_l += rho;
	  rb1_l += rho * Sigma_ff[i][i] / 2.0;
	  rb2_l += rho * PI * Sigma_ff[i][i] * Sigma_ff[i][i];
	  rb3_l += rho * (PI / 6.0) * Sigma_ff[i][i]
	    * Sigma_ff[i][i] * Sigma_ff[i][i];

	  if (Lsteady_state == TRUE)
	    rho = Rho_b_RTF[i];
	  else
	    rho = Rho_coex[0];

	  rb0_r += rho;
	  rb1_r += rho * Sigma_ff[i][i] / 2.0;
	  rb2_r += rho * PI * Sigma_ff[i][i] * Sigma_ff[i][i];
	  rb3_r += rho * (PI / 6.0) * Sigma_ff[i][i]
	    * Sigma_ff[i][i] * Sigma_ff[i][i];
	}
    }

  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3 * inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq * inv_one_m_rb3;

  dphi_drb_bulk_ptr->S0 = log (inv_one_m_rb3);
  dphi_drb_bulk_ptr->S1 = rb2 * inv_one_m_rb3;
  dphi_drb_bulk_ptr->S2 = rb1 * inv_one_m_rb3 +
    rb2 * rb2 * inv_one_m_rb3_sq / (8.0 * PI);
  dphi_drb_bulk_ptr->S3 = rb0 * inv_one_m_rb3 +
    rb1 * rb2 * inv_one_m_rb3_sq +
    rb2 * rb2 * rb2 * inv_one_m_rb3_3rd / (12.0 * PI);

  /* dphi_drb_bulk_ptr->V1 = dphi_drb_bulk_ptr->V2 = 0 in the bulk */
  /* -- not calculated */

  if (Iliq_vap == 3 || Lsteady_state)
    {
      inv_one_m_rb3 = 1.0 / (1.0 - rb3_l);
      inv_one_m_rb3_sq = inv_one_m_rb3 * inv_one_m_rb3;
      inv_one_m_rb3_3rd = inv_one_m_rb3_sq * inv_one_m_rb3;

      dphi_drb_bulk_left_ptr->S0 = log (inv_one_m_rb3);
      dphi_drb_bulk_left_ptr->S1 = rb2_l * inv_one_m_rb3;
      dphi_drb_bulk_left_ptr->S2 = rb1_l * inv_one_m_rb3 +
	rb2_l * rb2_l * inv_one_m_rb3_sq / (8.0 * PI);
      dphi_drb_bulk_left_ptr->S3 = rb0_l * inv_one_m_rb3 +
	rb1_l * rb2_l * inv_one_m_rb3_sq +
	rb2_l * rb2_l * rb2_l * inv_one_m_rb3_3rd / (12.0 * PI);

      inv_one_m_rb3 = 1.0 / (1.0 - rb3_r);
      inv_one_m_rb3_sq = inv_one_m_rb3 * inv_one_m_rb3;
      inv_one_m_rb3_3rd = inv_one_m_rb3_sq * inv_one_m_rb3;

      dphi_drb_bulk_right_ptr->S0 = log (inv_one_m_rb3);
      dphi_drb_bulk_right_ptr->S1 = rb2_r * inv_one_m_rb3;
      dphi_drb_bulk_right_ptr->S2 = rb1_r * inv_one_m_rb3 +
	rb2_r * rb2_r * inv_one_m_rb3_sq / (8.0 * PI);
      dphi_drb_bulk_right_ptr->S3 = rb0_r * inv_one_m_rb3 +
	rb1_r * rb2_r * inv_one_m_rb3_sq +
	rb2_r * rb2_r * rb2_r * inv_one_m_rb3_3rd / (12.0 * PI);
    }
}

/****************************************************************************/
/* d2phi_drb2_delta:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                 for the dphi_drb that use Delta_Fn Stencils (all but S3) */

struct RB_Struct
d2phi_drb2_delta (struct RB_Struct rho_bar, double weight,
		  int *offset, double inv_rad, double inv_4pir,
		  double inv_4pirsq, double *sign, int fill_flag)
{
  struct RB_Struct tmp;
  double rb0, rb1, rb2, rb3, rb1v[NDIM_MAX], rb2v[NDIM_MAX];
  double inv_one_m_rb3, inv_one_m_rb3_sq, inv_one_m_rb3_3rd;
  double vector[NDIM_MAX];
  int idim;

  rb3 = rho_bar.S3;
  rb2 = rho_bar.S2;
  rb1 = rho_bar.S1;
  rb0 = rho_bar.S0;

  for (idim = 0; idim < Ndim; idim++)
    {
      rb2v[idim] = sign[idim] * rho_bar.V2[idim];
      rb1v[idim] = sign[idim] * rho_bar.V1[idim];
    }

  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3 * inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq * inv_one_m_rb3;
  for (idim = 0; idim < Ndim; idim++)
    vector[idim] = offset[idim] * Esize_x[idim] * inv_rad;

  tmp.S2 = (inv_one_m_rb3 * inv_4pir
	    + rb2 * Inv_4pi * inv_one_m_rb3_sq) * weight;
  for (idim = 0; idim < Ndim; idim++)
    tmp.S2 += weight * vector[idim] * rb2v[idim] * Inv_4pi * inv_one_m_rb3_sq;

  tmp.S3 = weight * (inv_4pirsq * inv_one_m_rb3
		     + rb2 * inv_4pir * inv_one_m_rb3_sq
		     + rb1 * inv_one_m_rb3_sq +
		     rb2 * rb2 * Inv_4pi * inv_one_m_rb3_3rd);
  for (idim = 0; idim < Ndim; idim++)
    tmp.S3 += weight * inv_one_m_rb3_sq *
      (-rb2v[idim] * rb2v[idim] * Inv_4pi * inv_one_m_rb3
       + vector[idim]
       * (rb2v[idim] * inv_4pir + rb1v[idim]
	  + rb2 * rb2v[idim] * inv_one_m_rb3 / (2 * PI)));

  if (fill_flag != JAC_SAVE_FILL_1)
    {

      tmp.S0 = 0.0;
      tmp.S1 = inv_one_m_rb3 * weight;

      if (fill_flag == FULL_MSR_FILL)
	{

	  for (idim = 0; idim < Ndim; idim++)
	    tmp.V1[idim] =
	      sign[idim] * (weight * inv_one_m_rb3) * vector[idim];
	  for (idim = 0; idim < Ndim; idim++)
	    tmp.V2[idim] =
	      sign[idim] * weight * (-rb2v[idim] * Inv_4pi *
				     inv_one_m_rb3_sq -
				     (-inv_4pir * inv_one_m_rb3 -
				      rb2 * Inv_4pi * inv_one_m_rb3_sq) *
				     vector[idim]);
	}
    }

  return (tmp);
}

/****************************************************************************/
/* d2phi_drb2_delta2:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                 for the dphi_drb that use Delta_Fn Stencils (all but S3) */

struct RB_Struct
d2phi_drb2_delta2 (struct RB_Struct rho_bar, double weight,
		   int *offset, double inv_rad, double inv_4pir,
		   double inv_4pirsq, double *sign, int fill_flag)
{
  struct RB_Struct tmp;
  double rb0, rb1, rb2, rb3, rb1v[NDIM_MAX], rb2v[NDIM_MAX];
  double inv_one_m_rb3, inv_one_m_rb3_sq, inv_one_m_rb3_3rd;
  double vector[NDIM_MAX];
  int idim;
  double DOT_rho22, alpha, alpha_sq, alpha_cb, beta, gamma[3], DOT_gamma, eps;

  rb3 = rho_bar.S3;
  rb2 = rho_bar.S2;
  rb1 = rho_bar.S1;
  rb0 = rho_bar.S0;

  for (idim = 0; idim < Ndim; idim++)
    {
      rb2v[idim] = sign[idim] * rho_bar.V2[idim];
      rb1v[idim] = sign[idim] * rho_bar.V1[idim];
    }

  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3 * inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq * inv_one_m_rb3;
  for (idim = 0; idim < Ndim; idim++)
    vector[idim] = offset[idim] * Esize_x[idim] * inv_rad;

  DOT_rho22 = 0.0;
  for (idim = 0; idim < Ndim; idim++)
    {
      DOT_rho22 += rb2v[idim] * rb2v[idim];
    }
  if (rb2 > 1.e-15)
    {
      alpha = rb2 - DOT_rho22 / rb2;
      beta = 1.0 + DOT_rho22 / (rb2 * rb2);
      DOT_gamma = 0.0;
      for (idim = 0; idim < Ndim; idim++)
	{
	  gamma[idim] = rb2v[idim] / rb2;
	  DOT_gamma += gamma[idim] * gamma[idim];
	}
      eps = alpha / rb2;
    }
  else
    {
      alpha = rb2;
      beta = 1.0;
      for (idim = 0; idim < Ndim; idim++)
	gamma[idim] = 0.0;
      DOT_gamma = 0.0;
      eps = 1.0;
    }
  alpha_sq = alpha * alpha;
  alpha_cb = alpha_sq * alpha;

  tmp.S2 = inv_one_m_rb3 * inv_4pir * weight;	/*old rosen rb2-rb1 */

  tmp.S2 += weight * Inv_4pi * inv_one_m_rb3_sq *	/* new rosen rb2-rb2 */
    alpha * (beta * beta - eps * DOT_gamma);

  for (idim = 0; idim < Ndim; idim++)	/*new rosen rb2-rb2v */
    tmp.S2 -= weight * Inv_4pi * inv_one_m_rb3_sq *
      alpha * (eps - 2 * beta) * gamma[idim] * vector[idim];

  tmp.S3 = weight * (inv_4pirsq * inv_one_m_rb3	/* old rosen */
		     + rb2 * inv_4pir * inv_one_m_rb3_sq
		     + rb1 * inv_one_m_rb3_sq);

  for (idim = 0; idim < Ndim; idim++)	/* old rosen rb3-rbv1/v2 */
    tmp.S3 += weight * inv_one_m_rb3_sq *
      vector[idim] * (rb2v[idim] * inv_4pir + rb1v[idim]);

  tmp.S3 += weight * Inv_4pi * inv_one_m_rb3_3rd * alpha_sq * beta;	/*new rosen rb3-rb2 */

  for (idim = 0; idim < Ndim; idim++)
    tmp.S3 += weight * 2.0 * Inv_4pi * inv_one_m_rb3_3rd *	/*new rosen rb3-rb2v */
      alpha_sq * gamma[idim] * vector[idim];

  if (fill_flag != JAC_SAVE_FILL_1)
    {
      tmp.S0 = 0.0;		/*old rosen rb0-(rb0-rb2,rbv1,rbv2) */
      tmp.S1 = inv_one_m_rb3 * weight;	/*old rosen rb1-rb2 */

      if (fill_flag == FULL_MSR_FILL)
	{

	  for (idim = 0; idim < Ndim; idim++)
	    {			/*old rosen */
	      tmp.V1[idim] = sign[idim] * (weight * inv_one_m_rb3) * vector[idim];	/*rb1v-rb2v */
	      tmp.V2[idim] = sign[idim] * (weight * inv_4pir * inv_one_m_rb3) * vector[idim];	/*rb2v-rb1v */
	    }

	  for (idim = 0; idim < Ndim; idim++)
	    {
	      tmp.V2[idim] += sign[idim] * weight * Inv_4pi * inv_one_m_rb3_sq *	/*new rb2v-rb2 */
		alpha * (eps - 2.0 * beta) * gamma[idim];

	      tmp.V2[idim] -= sign[idim] * weight * Inv_4pi * inv_one_m_rb3_sq *	/*new rb2v-rb2v */
		alpha * (4.0 * DOT_gamma - eps) * vector[idim];
	    }
	}
    }

  return (tmp);
}

/****************************************************************************/
/* d2phi_drb2_theta:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                    for the dphi_drb that use Theta_Fn Stencils (S3)      */

struct RB_Struct
d2phi_drb2_theta (struct RB_Struct rho_bar, double weight,
		  int *offset, int fill_flag)
{
  struct RB_Struct tmp;
  double rb0, rb1, rb2, rb3, rb1v[NDIM_MAX], rb2v[NDIM_MAX];
  double inv_one_m_rb3, inv_one_m_rb3_sq, inv_one_m_rb3_3rd,
    inv_one_m_rb3_4th;
  int idim;

  rb3 = rho_bar.S3;
  rb2 = rho_bar.S2;
  rb1 = rho_bar.S1;
  rb0 = rho_bar.S0;

  for (idim = 0; idim < Ndim; idim++)
    {
      rb2v[idim] = rho_bar.V2[idim];
      rb1v[idim] = rho_bar.V1[idim];
    }
  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3 * inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq * inv_one_m_rb3;
  inv_one_m_rb3_4th = inv_one_m_rb3_3rd * inv_one_m_rb3;

  tmp.S2 = weight * (rb1 * inv_one_m_rb3_sq
		     + rb2 * rb2 * Inv_4pi * inv_one_m_rb3_3rd);
  for (idim = 0; idim < Ndim; idim++)
    tmp.S2 += -weight * rb2v[idim] * rb2v[idim] * Inv_4pi * inv_one_m_rb3_3rd;
  tmp.S3 = weight * (rb0 * inv_one_m_rb3_sq
		     + 2 * rb1 * rb2 * inv_one_m_rb3_3rd
		     + rb2 * rb2 * rb2 * Inv_4pi * inv_one_m_rb3_4th);
  for (idim = 0; idim < Ndim; idim++)
    tmp.S3 += -weight * (2 * rb1v[idim] * rb2v[idim] * inv_one_m_rb3_3rd
			 +
			 3 * rb2 * rb2v[idim] * rb2v[idim] * Inv_4pi *
			 inv_one_m_rb3_4th);

  if (fill_flag != JAC_SAVE_FILL_1)
    {
      tmp.S0 = weight * inv_one_m_rb3;
      tmp.S1 = weight * rb2 * inv_one_m_rb3_sq;

      if (fill_flag == FULL_MSR_FILL)
	{
	  for (idim = 0; idim < Ndim; idim++)
	    tmp.V1[idim] = -weight * rb2v[idim] * inv_one_m_rb3_sq;
	  for (idim = 0; idim < Ndim; idim++)
	    tmp.V2[idim] = -weight * (rb1v[idim] * inv_one_m_rb3_sq
				      +
				      rb2 * rb2v[idim] * inv_one_m_rb3_3rd /
				      (2 * PI));
	}
    }

  return (tmp);
}

/****************************************************************************/
/* d2phi_drb2_theta2:  calculate the derivatives of the dphi_drb w.r.t. rb   */
/*                    for the dphi_drb that use Theta_Fn Stencils (S3)      */

struct RB_Struct
d2phi_drb2_theta2 (struct RB_Struct rho_bar, double weight,
		   int *offset, int fill_flag)
{
  struct RB_Struct tmp;
  double rb0, rb1, rb2, rb3, rb1v[NDIM_MAX], rb2v[NDIM_MAX];
  double inv_one_m_rb3, inv_one_m_rb3_sq, inv_one_m_rb3_3rd,
    inv_one_m_rb3_4th;
  int idim;
  double DOT_rho22, DOT_rho12, alpha, alpha_sq, alpha_cb, beta, gamma[3];

  rb3 = rho_bar.S3;
  rb2 = rho_bar.S2;
  rb1 = rho_bar.S1;
  rb0 = rho_bar.S0;

  for (idim = 0; idim < Ndim; idim++)
    {
      rb2v[idim] = rho_bar.V2[idim];
      rb1v[idim] = rho_bar.V1[idim];
    }
  inv_one_m_rb3 = 1.0 / (1.0 - rb3);
  inv_one_m_rb3_sq = inv_one_m_rb3 * inv_one_m_rb3;
  inv_one_m_rb3_3rd = inv_one_m_rb3_sq * inv_one_m_rb3;
  inv_one_m_rb3_4th = inv_one_m_rb3_3rd * inv_one_m_rb3;

  DOT_rho22 = 0.0;
  DOT_rho12 = 0.0;
  for (idim = 0; idim < Ndim; idim++)
    {
      DOT_rho22 += rb2v[idim] * rb2v[idim];
      DOT_rho12 += rb1v[idim] * rb2v[idim];
    }

  if (rb2 > 1.e-15)
    {
      alpha = rb2 - DOT_rho22 / rb2;
      beta = 1.0 + DOT_rho22 / (rb2 * rb2);
      for (idim = 0; idim < Ndim; idim++)
	{
	  gamma[idim] = rb2v[idim] / rb2;
	}
    }
  else
    {
      alpha = rb2;
      beta = 1.0;
      for (idim = 0; idim < Ndim; idim++)
	gamma[idim] = 0.0;
    }
  alpha_sq = alpha * alpha;
  alpha_cb = alpha_sq * alpha;

  tmp.S2 = weight * rb1 * inv_one_m_rb3_sq;	/*old rb2-rb3 */
  tmp.S2 += weight * Inv_4pi * inv_one_m_rb3_3rd * alpha_sq * beta;	/*new rb2-rb
									   3 */

  tmp.S3 = weight * (rb0 * inv_one_m_rb3_sq	/*old rb3-rb3 */
		     + 2.0 * (rb1 * rb2 - DOT_rho12) * inv_one_m_rb3_3rd);
  tmp.S3 += weight * Inv_4pi * inv_one_m_rb3_4th * alpha_cb;	/*new rb3-rb3 */

  if (fill_flag != JAC_SAVE_FILL_1)
    {
      tmp.S0 = weight * inv_one_m_rb3;	/*old rb0-rb3 */
      tmp.S1 = weight * rb2 * inv_one_m_rb3_sq;	/*old rb1-rb3 */

      if (fill_flag == FULL_MSR_FILL)
	{
	  for (idim = 0; idim < Ndim; idim++)
	    {
	      tmp.V1[idim] = -weight * rb2v[idim] * inv_one_m_rb3_sq;	/*old rbv1-rb3 */
	      tmp.V2[idim] = -weight * rb1v[idim] * inv_one_m_rb3_sq;	/*old rbv2-rb3 */
	    }

	  for (idim = 0; idim < Ndim; idim++)	/*new rbv2-rb3 */
	    tmp.V2[idim] -=
	      weight * 2.0 * Inv_4pi * inv_one_m_rb3_3rd * alpha_sq *
	      gamma[idim];

	}
    }

  return (tmp);
}

/****************************************************************************/
/* load_Jac_Save_Jacobian:  assemble the Hard Spheres  jacobian !!!                  */

void
load_Jac_Save_Jacobian (double *mat_row, struct RB_Struct tmp,
			int jnode_box, int jzone, int *bindx_tmp,
			int fill_flag, double ***jac_weights_hs,
			int ***jac_columns_hs)
{
  int j, jn_loc, *jac_col, isten, njcol;
  double *jac_wt, d2;

  jn_loc = B2L_1stencil[jnode_box];

  for (isten = 0; isten <= 1; isten++)
    {

      jac_col = jac_columns_hs[jn_loc][isten];
      njcol = jac_col[0];

      if (fill_flag != MSR_PREPROCESS)
	{

	  jac_wt = jac_weights_hs[jn_loc][isten];

	  if (isten == 0)
	    {
	      d2 = tmp.S2;
	      if (fill_flag == JAC_SAVE_FILL_2)
		d2 += tmp.S0 * Inv_4pirsq[0] + tmp.S1 * Inv_4pir[0];
	    }
	  else
	    d2 = tmp.S3;

	  /* regular fill, jac_col holds local unknown numbers */

	  for (j = 1; j < njcol; j++)
	    {
	      mat_row[*(++jac_col)] += d2 * *(++jac_wt);
	    }
	}
      else
	{

	  /* preprocess fill, jac_col holds global unknown numbers */

	  for (j = 1; j < njcol; j++)
	    bindx_tmp[*(++jac_col)] = TRUE;
	}
    }
  return;
}

/****************************************************************************/
/* load_Jacobian:  assemble the Hard Spheres  jacobian !!!                  */

void
load_Jacobian (double *mat_row, struct RB_Struct tmp, int jnode_box,
	       int jzone, int *bindx_tmp, int fast_fill_TF_j, int fill_flag)
{
  int idim, jlist = 0;
  int jcomp;
  double tmp_scalar;
  double tmp_vector[3];

  /*
   * Loop over Delta_Fn Stencil to assemble the following rho_bars: 
   * S0, S1, S2, V1[Ndim], V2[Ndim]
   */

  for (jcomp = 0; jcomp < Ncomp; jcomp++)
    {
      if (Nlists_HW == 1 || Nlists_HW == 2)
	jlist = 0;
      else if (Nlists_HW > 2)
	jlist = jcomp;

      /* precalculate quantities used in fill below */

      if (fill_flag != MSR_PREPROCESS)
	{

	  tmp_scalar = tmp.S0 * Inv_4pirsq[jcomp]
	    + tmp.S1 * Inv_4pir[jcomp] + tmp.S2;

	  for (idim = 0; idim < Ndim; idim++)
	    tmp_vector[idim] = Esize_x[idim] * Inv_rad[jcomp] *
	      (tmp.V1[idim] * Inv_4pir[jcomp] + tmp.V2[idim]);

	}
      else
	{
	  tmp_scalar = 0.0;
	  tmp_vector[0] = 0.0;
	  tmp_vector[1] = 0.0;
	  tmp_vector[2] = 0.0;
	}

      /* Pick which fill depending on Fast_fill_TF flag */

      switch (fast_fill_TF_j)
	{
	case CHECK_NONE:
	  fast_fill_Jacobian (DELTA_FN, jcomp, jzone, jnode_box,
			      tmp_scalar, tmp_vector,
			      mat_row, bindx_tmp, fill_flag);
	  break;
	case CHECK_HW:
	  medium_fill_Jacobian (DELTA_FN, jcomp, jzone, jlist, jnode_box,
				tmp_scalar, tmp_vector,
				mat_row, bindx_tmp, fill_flag);
	  break;
	case CHECK_BC:
	case CHECK_BOTH:
	  slow_fill_Jacobian (DELTA_FN, jcomp, jzone,
			      jlist, jnode_box, tmp_scalar, tmp_vector,
			      mat_row, bindx_tmp, fill_flag);
	  break;
	}

      /* Now do the same for the lone rho_bar that needs the Theta_Fn  */

      if (fill_flag != MSR_PREPROCESS)
	tmp_scalar = tmp.S3;

      /* Pick which fill depending on Fast_fill_TF flag */

      switch (fast_fill_TF_j)
	{
	case CHECK_NONE:
	  fast_fill_Jacobian (THETA_FN, jcomp, jzone, jnode_box,
			      tmp_scalar, tmp_vector,
			      mat_row, bindx_tmp, fill_flag);
	  break;
	case CHECK_HW:
	  medium_fill_Jacobian (THETA_FN, jcomp, jzone, jlist, jnode_box,
				tmp_scalar, tmp_vector,
				mat_row, bindx_tmp, fill_flag);
	  break;
	case CHECK_BC:
	case CHECK_BOTH:
	  slow_fill_Jacobian (THETA_FN, jcomp, jzone,
			      jlist, jnode_box, tmp_scalar, tmp_vector,
			      mat_row, bindx_tmp, fill_flag);
	  break;
	}
    }
}

/****************************************************************************/

#define FINAL_JACOBIAN_FILL_MACRO                                            \
{                                                                            \
  if (fill_flag != MSR_PREPROCESS){                                          \
     loc_k = B2L_unknowns[knode_box*Nunk_per_node + jcomp];                  \
     if (sten_type == DELTA_FN){                                             \
        switch (Ndim) {                                                      \
         case 1:  mat_row[loc_k] +=                                          \
                    weight * (tmp_scalar + offset[0] * tmp_vector[0]); break;\
                                                                             \
         case 2:  mat_row[loc_k] +=                                          \
                    weight * (tmp_scalar + offset[0] * tmp_vector[0] +       \
                                           offset[1] * tmp_vector[1]); break;\
         case 3:  mat_row[loc_k] +=                                          \
                    weight * (tmp_scalar + offset[0] * tmp_vector[0] +       \
                                           offset[1] * tmp_vector[1] +       \
                                           offset[2] * tmp_vector[2]);       \
        }                                                                    \
      }                                                                      \
      else{ /*THETA_FN*/                                                     \
                  mat_row[loc_k] += tmp_scalar * weight;                     \
      }                                                                      \
   }                                                                         \
   else{  /* fill_flag == MSR_PREPROCESS */                                  \
      bindx_tmp[knode_box*Nunk_per_node + jcomp] = TRUE;                     \
   }                                                                         \
}				/* end of fast_Jacobian_fill macro */

/****************************************************************************/
/*fast_fill_Jacobian: This routine contains the logic for a fast fill where
                     we don't have to worry about hard walls or boundaries. */

void
fast_fill_Jacobian (int sten_type,
		    int jcomp, int jzone, int jnode_box,
		    double tmp_scalar, double *tmp_vector,
		    double *mat_row, int *bindx_tmp, int fill_flag)
{
  double *sten_weight, weight;
  int **sten_offset, *offset;
  int knode_box, isten, loc_k;
  struct Stencil_Struct *sten;

  sten = &(Stencil[sten_type][jzone][jcomp]);
  sten_offset = sten->Offset;
  sten_weight = sten->Weight;

  for (isten = 0; isten < sten->Length; isten++)
    {

      /* get offset and weight for this stencil point */

      offset = sten_offset[isten];
      weight = sten_weight[isten];

      /* Instead of offset_to_node, do this quick calc */

      knode_box = jnode_box;
      switch (Ndim)
	{
	case 3:
	  knode_box += offset[0] + offset[1] * Nodes_x_box[0]
	    + offset[2] * Nodes_plane_box;
	  break;
	case 2:
	  knode_box += offset[0] + offset[1] * Nodes_x_box[0];
	  break;
	case 1:
	  knode_box += offset[0];
	}

      if (!Zero_density_TF[knode_box][jcomp])
	{

	FINAL_JACOBIAN_FILL_MACRO}

    }				/* End of loop over Stencil Length */

  return;
}

/****************************************************************************/
/*medium_fill_Jacobian: This routine contains the logic for a medium fill
                           where worry about hard walls but not boundaries. */

void
medium_fill_Jacobian (int sten_type,
		      int jcomp, int jzone, int jlist, int jnode_box,
		      double tmp_scalar, double *tmp_vector,
		      double *mat_row, int *bindx_tmp, int fill_flag)
{
  double *sten_weight, weight;
  int **sten_offset, *offset;
  int knode_box, isten, loc_k;
  struct Stencil_Struct *sten;
  int reflect_flag[3] = { FALSE, FALSE, FALSE };

  sten = &(Stencil[sten_type][jzone][jcomp]);
  sten_offset = sten->Offset;
  sten_weight = sten->Weight;

  for (isten = 0; isten < sten->Length; isten++)
    {

      /* get offset and weight for this stencil point */

      offset = sten_offset[isten];
      weight = sten_weight[isten];

      /* Instead of offset_to_node, do this quick calc */

      knode_box = jnode_box;
      switch (Ndim)
	{
	case 3:
	  knode_box += offset[0] + offset[1] * Nodes_x_box[0]
	    + offset[2] * Nodes_plane_box;
	  break;
	case 2:
	  knode_box += offset[0] + offset[1] * Nodes_x_box[0];
	  break;
	case 1:
	  knode_box += offset[0];
	}

      if (!Zero_density_TF[knode_box][jcomp])
	{

	  if (fill_flag != MSR_PREPROCESS)
	    if (Nodes_2_boundary_wall[jlist][knode_box] != -1)
	      weight = HW_boundary_weight
		(jcomp, jlist, sten->HW_Weight[isten], knode_box,
		 reflect_flag);


	FINAL_JACOBIAN_FILL_MACRO}

    }				/* End of loop over Stencil Length */

  return;
}

/****************************************************************************/
/*slow_fill_Jacobian: This routine contains the logic for a slow fill where
                      we worry about both hard walls and boundaries.        */

void
slow_fill_Jacobian (int sten_type,
		    int jcomp, int jzone, int jlist, int jnode_box,
		    double tmp_scalar, double *tmp_vector,
		    double *mat_row, int *bindx_tmp, int fill_flag)
{
  double *sten_weight, weight;
  int **sten_offset, *offset;
  int knode_box, isten, ijk_jnode_box[3], loc_k;
  int reflect_flag[NDIM_MAX];
  struct Stencil_Struct *sten;

  node_box_to_ijk_box (jnode_box, ijk_jnode_box);

  sten = &(Stencil[sten_type][jzone][jcomp]);
  sten_offset = sten->Offset;
  sten_weight = sten->Weight;

  for (isten = 0; isten < sten->Length; isten++)
    {

      /* get offset and weight for this stencil point */

      offset = sten_offset[isten];
      weight = sten_weight[isten];

      /* Find in the Stencil position on mesh */

      knode_box = offset_to_node_box (ijk_jnode_box, offset, reflect_flag);

      if (knode_box > -1)
	{

	  if (!Zero_density_TF[knode_box][jcomp])
	    {

	      /*
	       * If we are doing Hard Walls, and this is a boundary node, 
	       * modify the weight to only include contributions from fluid elements
	       * and not wall elements (where the density is 0 but, through linear
	       * interpolation, would otherwise give a nonzero contribution).
	       */

	      if (Lhard_surf && fill_flag != MSR_PREPROCESS)
		if (Nodes_2_boundary_wall[jlist][knode_box] != -1)
		  weight = HW_boundary_weight
		    (jcomp, jlist, sten->HW_Weight[isten], knode_box,
		     reflect_flag);

	    FINAL_JACOBIAN_FILL_MACRO}
	}

    }				/* End of loop over Stencil Length */

  return;
}
