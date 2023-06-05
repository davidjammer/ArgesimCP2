#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
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
	exit(code);
}

int main(int argc, char *argv[])
{
	double **geometry, **rho, **ux, **uy, **usqr, **u;
	double ***f, ***feq;
	double viscosity, tau;
	index2d *driving_ids; 
	int n_driving_ids;
	int rep;
	int x, y;
	char filename[128];
	
	int Ny = DEFAULT_Ny;
	int Nx = DEFAULT_Ny;
	int REP = DEFAULT_REP;
	
	if(set_nx_from_argv(argc, argv, &Nx) < 0)
	{
		
		printf("Incorrect number for Nx\n");
		lbm_exit(EXIT_FAILURE);
	}
	if(set_ny_from_argv(argc, argv, &Ny) < 0)
	{
		
		printf("Incorrect number for Ny\n");
		lbm_exit(EXIT_FAILURE);
	}
	if(set_rep_from_argv(argc, argv, &REP) < 0)
	{
		
		printf("Incorrect number for Rep\n");
		lbm_exit(EXIT_FAILURE);
	}
	
	sprintf(filename, "u.dat");
	set_filename_from_argv(argc, argv, filename);
	
	if(print_help(argc, argv) == 1)
		lbm_exit(EXIT_SUCCESS);

	//print parameter
	
	printf("Nx: %d\n", Nx);
	printf("Ny: %d\n", Ny);
	printf("Rep: %d\n", REP);
	printf("Output file: %s\n", filename);
	

	
	//calculate the viscosity [Equation 6]
	viscosity = (Ny - 1) * u_0 / RE;
	//calculate tau [Equation 5]
	tau = (6 * viscosity + 1) / 2;
	
	//create the geometry
	geometry = create_matrix2d(Nx, Ny);	
	create_geometry(geometry, Nx, Ny);
	get_driving_ids(geometry, Nx, Ny, &driving_ids, &n_driving_ids);

	//create the 2d matrix for rho, ux, uy and usqr
	rho = create_matrix2d(Nx, Ny);
	ux = create_matrix2d(Nx, Ny);
	uy = create_matrix2d(Nx, Ny);
	usqr = create_matrix2d(Nx, Ny);
	
	//create the 3d matrix for fand feq
	f = create_matrix3d(Nx, Ny, 9);
	feq = create_matrix3d(Nx, Ny, 9);
	
	//initialize the distribution values
	init_dist_func_values(f, Nx, Ny);
	
	//simualtion loop
	for(rep=0; rep<REP; rep++)
	{
		//collision step
		calc_rho(rho, f, Nx, Ny);
		calc_ux(ux, f, rho, Nx, Ny);
		calc_uy(uy, f, rho, Nx, Ny);
		
		set_ux_for_driving_cells(ux, driving_ids, n_driving_ids);
		set_uy_for_driving_cells(uy, driving_ids, n_driving_ids);
		calc_usqr(usqr, ux, uy, Nx, Ny);
		
		clac_feq(feq, rho, ux, uy, usqr, Nx, Ny); 
		
		f_cell_calc(f, feq, geometry, tau, Nx, Ny);
		
		//propagation step
		propagation_step(f, Nx, Ny);

#ifdef DEBUG
		printf("Iteration %d/%d\n", rep, REP);
#endif		
	}
	
	//calculate the relative macroscopic velocity magnitude and save it
	u = calc_magnitude(Nx, Ny, ux, uy);
	
	save(filename, u, Nx, Ny);

	return EXIT_SUCCESS;
}
