/*
  Prototypes from dft_out_energy.c
*/
double calc_free_energy (FILE * fp, double *x, double fac_area,
			 double fac_vol, int print_flag);

void assemble_HS_free_energy (double *x, double *sum_phispt,
			      double *sum_phispt_b, double *sum_phispt_b_old);

double energy_elec (double *x, double *sum3);

double charge_stress (double *x, double *sum2);

double free_energy_charging_up (double *x);

double int_stencil (double *x, int inode_box, int icomp, int sten_type);

double phispt_i (struct RB_Struct *rho_bar, int inode);

double phispt_bulk ();

double calc_deriv_e (int idim, int inode0, int flag, int *blocked, double *x,
		     int ilist);

double calc_deriv2 (int idim, int inode0, int flag, double *x);

double calc_free_energy_polymer (FILE * fp, double *x, double fac_area,
				 double fac_vol);

double calc_u_ideal (int itype_mer, int *ijk_box, double *x, double *fluid,
		     double *bulk);
