/*
  Prototypes for dft_mesh_surfaces.c
*/

void setup_surface (FILE * fp2, int *nelems_f,
		    int **nelems_w_per_w,
		    int **elems_f,
		    int ***elems_w_per_w, int *elem_zones,
		    int *fast_fill_elem_TF, int ***el_type);

void find_wall_images (int idim, int *image, double **image_pos, double *pos);

void flag_wall_el (int inode, int ilist, int iwall, int iel_box, int **L_wall,
		   int **nelems_w_per_w, int ***elems_w_per_w,
		   int ***el_type);

void els_bumpy (int iwall, int real_wall, int itype, int **L_wall,
		double *x_min, int **nelems_w_per_w, int ***elems_w_per_w,
		FILE * fp2, int *fast_fill_elem_TF, int ***el_type,
		double **image_pos);

void els_spheres (int iwall, int real_wall, int itype, int **L_wall,
		  double *x_min, int **nelems_w_per_w, int ***elems_w_per_w,
		  int *fast_fill_elem_TF, int ***el_type, double **image_pos);

void els_cyls_3D (int iwall, int real_wall, int itype, int **L_wall,
		  double *x_min, int **nelems_w_per_w, int ***elems_w_per_w,
		  int *fast_fill_elem_TF, int ***el_type, double **image_pos);

void els_cyls_cos_3D (int iwall, int real_wall, int itype, int **L_wall,
		      double *x_min, int **nelems_w_per_w,
		      int ***elems_w_per_w, int *fast_fill_elem_TF,
		      int ***el_type, double **image_pos);

void els_atomic_centers (int iwall, int real_wall, int itype, int **L_wall,
			 double *x_min, int **nelems_w_per_w,
			 int ***elems_w_per_w, int *fast_fill_elem_TF,
			 int ***el_type, double **image_pos);

void els_slit_pore_2D (int iwall, int real_wall, int itype, int **L_wall,
		       double *x_min, int **nelems_w_per_w,
		       int ***elems_w_per_w, int *fast_fill_elem_TF,
		       int ***el_type, double **image_pos);

void els_cyl_pore_3D (int iwall, int real_wall, int itype, int **L_wall,
		      double *x_min, int **nelems_w_per_w,
		      int ***elems_w_per_w, int *fast_fill_elem_TF,
		      int ***el_type, double **image_pos);

void els_cone_pore_2D (int iwall, int real_wall, int itype, int **L_wall,
		       double *x_min, int **nelems_w_per_w,
		       int ***elems_w_per_w, int *fast_fill_elem_TF,
		       int ***el_type, double **image_pos);

void els_cone_pore_3D (int iwall, int real_wall, int itype, int **L_wall,
		       double *x_min, int **nelems_w_per_w,
		       int ***elems_w_per_w, int *fast_fill_elem_TF,
		       int ***el_type, double **image_pos);

void els_cyl_pores (int iwall, int real_wall, int itype, int **L_wall,
		    double *x_min, int **nelems_w_per_w, int ***elems_w_per_w,
		    int *fast_fill_elem_TF, int ***el_type,
		    double **image_pos);

void els_planar (int iwall, int real_wall, int itype,
		 int **L_wall, double *x_min,
		 int **nelems_w_per_w, int ***elems_w_per_w,
		 int *fast_fill_elem_TF, int ***el_type, double **image_pos);

void els_finite_planar (int iwall, int real_wall, int itype,
			int **L_wall, double *x_min,
			int **nelems_w_per_w, int ***elems_w_per_w,
			int *fast_fill_elem_TF, int ***el_type,
			double **image_pos);
