/* 
   Prototypes for dft_out_main.c
*/

void post_process (double *x_internal,char *output_file3,int *niters,
                   double *time_save, int loop1, int binodal_flag);
void setup_integrals();
void print_header(char *output_file3);
void print_cont_variable(int cont_type,FILE *fp);
void print_cont_type(int cont_type,FILE *fp);
