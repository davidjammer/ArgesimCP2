#ifndef __LBM_MAGNITUDE
	#define __LBM_MAGNITUDE
	
	double **calc_magnitude(int Nx, int Ny, int nx, int ny, double **ux, double **uy, int dim_size[2], MPI_Comm lbm_comm);
	
#endif
