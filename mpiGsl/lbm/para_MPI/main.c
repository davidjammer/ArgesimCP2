#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

#include <lbm_parameter.h>
#include <functions.h>
#include <lbm_init.h>
#include <lbm_collision_step.h>
#include <lbm_propagation_step.h>
#include <lbm_magnitude.h>

//#define DEBUG

void lbm_exit(int code)
{
	MPI_Finalize();
	exit(code);
}

int main(int argc, char *argv[])
{
	MPI_Comm lbm_comm;
	double viscosity, tau;
	double **geometry, **rho, **ux, **uy, **usqr, **u;
	double ***f, ***feq;
	int world_size, world_rank;
	int nx, ny;
	index2d *driving_ids; 
	int n_driving_ids;
	int rep, x, y;
	char filename[128];
	
	int dim_size[2] = {0,0};
	int periods[2] = {0,0};
	int Ny = DEFAULT_Ny;
	int Nx = DEFAULT_Ny;
	int REP = DEFAULT_REP;
	
	//Initialization of MPI
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	
	if(IsPowerOfTwo(world_size) == 0)
	{
		if(world_rank == 0)
			printf("Incorrect number for processes. Number of processes musst be power of two.\n");
		lbm_exit(EXIT_FAILURE);
	}
	
	if(set_nx_from_argv(argc, argv, &Nx) < 0)
	{
		if(world_rank == 0)
			printf("Incorrect number for Nx\n");
		lbm_exit(EXIT_FAILURE);
	}
	if(set_ny_from_argv(argc, argv, &Ny) < 0)
	{
		if(world_rank == 0)
			printf("Incorrect number for Ny\n");
		lbm_exit(EXIT_FAILURE);
	}
	if(set_rep_from_argv(argc, argv, &REP) < 0)
	{
		if(world_rank == 0)
			printf("Incorrect number for Rep\n");
		lbm_exit(EXIT_FAILURE);
	}
	
	sprintf(filename, "u.dat");
	set_filename_from_argv(argc, argv, filename);
	
	if(print_help(argc, argv) == 1)
		lbm_exit(EXIT_SUCCESS);

	//print parameter
	if(world_rank == 0)
	{
		printf("Nx: %d\n", Nx);
		printf("Ny: %d\n", Ny);
		printf("Rep: %d\n", REP);
		printf("Output file: %s\n", filename);
	}
		
		
#ifdef DEBUG
	if(world_rank == 0)
		printf("world_size: %d\n", world_size);
#endif

	//calculate the viscosity [Equation 6]
	viscosity = (Ny - 1) * u_0 / RE;
	//calculate tau [Equation 5]
	tau = (6 * viscosity + 1) / 2;

	MPI_Dims_create(world_size, 2, dim_size);
	if((dim_size[0] * dim_size[1]) != world_size)
		lbm_exit(EXIT_FAILURE);
	
	//create communicator
	MPI_Cart_create(MPI_COMM_WORLD,2,dim_size,periods,0,&lbm_comm);
		
#ifdef DEBUG
	if(world_rank == 0)
		printf("dimension: [%d,%d]\n", dim_size[0], dim_size[1]);
#endif

	//get the size of local matrices
	ny = Ny / dim_size[0];
	nx = Nx / dim_size[1];

#ifdef DEBUG
	if(world_rank == 0)
		printf("local matrices size: [%d,%d]\n", nx, ny);
#endif

	//create the geometry
	geometry = dist_geometry(Nx, Ny, nx, ny, dim_size, lbm_comm);
	
	get_driving_ids(geometry, nx, ny, &driving_ids, &n_driving_ids);
	
	MPI_Barrier(lbm_comm);
	
	//create the 2d matrix for rho, ux, uy and usqr
	rho = create_matrix2d(nx, ny);
	ux = create_matrix2d(nx, ny);
	uy = create_matrix2d(nx, ny);
	usqr = create_matrix2d(nx, ny);

	//create the 3d matrix for f and feq
	f = create_matrix3d(nx, ny, 9);
	feq = create_matrix3d(nx, ny, 9);
	
	
	//initialize the distribution values
	init_dist_func_values(f, nx, ny);
	
	
	//simualtion loop
	for(rep=0; rep<REP; rep++)
	{
		//collision step
		calc_rho(rho, f, nx, ny);
		calc_ux(ux, f, rho, nx, ny);
		calc_uy(uy, f, rho, nx, ny);
		
		set_ux_for_driving_cells(ux, driving_ids, n_driving_ids);
		set_uy_for_driving_cells(uy, driving_ids, n_driving_ids);
		calc_usqr(usqr, ux, uy, nx, ny);
		
		calc_feq(feq, rho, ux, uy, usqr, nx, ny); 
		
		f_cell_calc(f, feq, geometry, tau, nx, ny);
		
		//propagation step
		propagation_step(f, nx, ny, lbm_comm);
		
#ifdef DEBUG
		if(world_rank == 0)
			printf("Iteration %d/%d\n", rep, REP);
#endif
	}
	
	//calculate the relative macroscopic velocity magnitude and save it
	u = calc_magnitude(Nx, Ny, nx, ny, ux, uy, dim_size, lbm_comm);
	
	if(u != NULL)
		save(filename, u, Nx, Ny);
	
	
	free_matrix2d(geometry);
	free_matrix2d(rho);
	free_matrix2d(ux);
	free_matrix2d(uy);
	free_matrix2d(usqr);
	if(u != NULL)
		free_matrix2d(u);
	
	free_matrix3d(f, nx, ny, 9);
	free_matrix3d(feq, nx, ny, 9);
	
	free(driving_ids);
	
	MPI_Finalize();
	
	return EXIT_SUCCESS;
}
