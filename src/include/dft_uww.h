/*
  Prototypes for dft_uww.c
*/
void setup_wall_wall_potentials( int **nelems_w_per_w,
                            int ***elems_w_per_w);
void setup_lj_atomic(int,int);
void setup_coulomb_atomic(int,int);
void setup_ww_integrated(int,int *,int **,int,int);
