/*
  Prototypes for dft_mesh.c
*/
void set_up_mesh (char *output_file1,char *output_file2);
void free_mesh_arrays(void);
void control_mesh(FILE *fp1,char *output_file2,int print_flag);
void setup_basic_domain(FILE *fp1);
void setup_basic_box(FILE *fp1);
void initialize_fast_fill_flag(void);
void boundary_setup(char *output_file1);
void boundary_free(void);
void setup_zeroTF_and_Node2bound (FILE *fp1,int ***el_type);
void setup_zeroTF_and_Node2bound_new (FILE *fp1,int ***el_type);
void boundary_properties(FILE *fp1);
void find_local_els(int inode,int *iel, int *iel_box,int flag);
void surf_el_to_list(int loc_inode, int ilist, int *iel_box,
		     int el, int type, int normal, int idim, 
		     double esize1,double esize2);
void setup_surface_charge(FILE *fp1);
void bc_setup_const_charge(int iwall_type, int loc_inode);
void setup_linear_grad_of_charge(void);
void setup_volume_charge1(int iwall);
void setup_volume_charge2(void);
void els_charge_spheres(double radius,double *x,int *nelems,
                        int *nelems_unique, int *elems,int charge_type);
void zones_el_to_nodes(int *elem_zones);
void fast_fill_el_to_nodes(int *fast_fill_el_TF);
void set_mesh_coarsen_flag(void);
void setup_area_IC(void);
void initialize_Aztec(void);
void MY_read_update(int *N_update, int *update[], int proc_config[],
                    int N, int *nodes_x, int chunk, int input_option);
