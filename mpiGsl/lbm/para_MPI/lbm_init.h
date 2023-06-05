
#ifndef __LB_INIT
	#define __LB_INIT
	
	double **dist_geometry(int Nx, int Ny, int nx, int ny, int dim_size[2], MPI_Comm lbm_comm);
	void create_geometry(double **g, int nx, int ny);
	void init_dist_func_values(double ***f, int nx, int ny);
#endif
