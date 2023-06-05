//mpicc -O3 main.c `pkg-config --libs --cflags gsl`
// call: para [REPS STEPS]

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <mpi.h>

#define REPS 1000   //number of simulations
#define STEPS 200   //number of compute steps

//parameter for the damped dynamic mass - spring system
#define T_END 2.0
#define K 9000.0
#define M 450.0
#define D_MIN 800.0
#define D_MAX 1200.0

//DGL for the damped dynamic mass - spring system
int f(double t, const double y[], double dydt[], void *params)
{
	double k, d, m;
	
	k = ((double*)params)[0];
	d = ((double*)params)[1];
	m = ((double*)params)[2];

	dydt[0] = y[1];
	dydt[1] =  - d / m * y[1] - k / m * y[0];
	
	return GSL_SUCCESS;
}

//save function; save t and y in file
void save(char *file, double *t, double *y, int n)
{
	FILE *fd;
	int i;
	
	fd=fopen(file,"w");
	
	for(i=0; i<n; i++)
	{
		fprintf(fd,"%E ",t[i]);
		fprintf(fd,"%E ",y[i]);
		
		fprintf(fd,"\n");
	}
	fclose(fd);

}

int main(int argc, char *argv[])
{
        long nReps;
	long nSteps;
        int world_size;
	int world_rank;
	int dim = 2;
	int rep, i;
	int l_rep;
	int step;
	double *tvec;
	double *y0vec, *y0vec_recv;
	const gsl_odeiv2_step_type *type_ptr;
	gsl_odeiv2_step *step_ptr;
	gsl_odeiv2_evolve *evolve_ptr;
	gsl_odeiv2_system my_system;
	const gsl_rng_type *T;
	gsl_rng *r;
	double y[2];
	double t;
	double tEnd;
	double params[3];

	// process command line
	if (argc <= 2)
	{
	  nReps = (long) REPS;
	  nSteps = (long) STEPS;
	} else {
	  nReps = atol(argv[1]);
	  nSteps = atol(argv[2]);
	}
	  
	//Initialisierung
	MPI_Init(NULL, NULL);
    	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	//vector for the time and the mean motion
	tvec = calloc(sizeof(double), nSteps+1);
	y0vec = calloc(sizeof(double), nSteps+1);

	type_ptr = gsl_odeiv2_step_rk4; // RK4 Solver
	step_ptr = gsl_odeiv2_step_alloc(type_ptr, dim); // stepping function
	evolve_ptr = gsl_odeiv2_evolve_alloc(dim); // evolution function
	T = gsl_rng_default; // default rng = mersenne twister
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, (unsigned long int) world_rank);
	        // different seed for each task
	
	//parameter vector for the DGL
	params[0] = K;
	params[2] = M; 
	tEnd = T_END;
	
	//DGL-System
	my_system.function = f;      // the right-hand-side functions dy[i]/dt
    	my_system.jacobian = NULL;   // the Jacobian df[i]/dy[j]
    	my_system.dimension = dim;   // number of diffeq's
    	my_system.params = params;   // parameters to pass to rhs and jacobian
	
	
	l_rep = nReps / world_size;
	if( (nReps % world_size) > 0 && (world_size - (nReps % world_size)) <= world_rank )
		l_rep++;
#ifdef DEBUG		
	printf("%3d: l_rep=%d\n", world_rank, l_rep);
#endif
	for(rep=0; rep<l_rep; rep++)
	{
		//set the state vector motion=0 and velocity=0.1 
		y[0] = 0.0;
		y[1] = 0.1;
		//set time
		t = 0.0;
	
		//set damping factor
		params[1] = gsl_ran_flat(r, D_MIN, D_MAX);
#ifdef DEBUG		
	        printf("%3d(%3d): d=%7.2f\n", rep, world_rank, params[1]);
#endif
		
		//first values
		tvec[0] = t;
		y0vec[0] += y[0];
		
		//sim. loop
		for(step=1; step <= nSteps; step++)
		{
			//ODE-Solver
			gsl_odeiv2_evolve_apply_fixed_step(evolve_ptr, NULL, step_ptr,
		                            &my_system, &t, tEnd/nSteps, y);

			//set result for time t
		    	tvec[step] = t;
			y0vec[step] += y[0];
		}
	
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(world_rank == 0)
	{
	        y0vec_recv = malloc(sizeof(double) * (nSteps+1));
		MPI_Reduce(y0vec, y0vec_recv, nSteps+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		//compute the mean motion
		for(step=0; step < nSteps+1; step++)
			y0vec_recv[step] /= nReps;

		//save results
		save("daten.dat", tvec, y0vec_recv, nSteps+1);
	}
	else
		MPI_Reduce(y0vec, NULL, nSteps+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
	//final steps
	gsl_odeiv2_evolve_free(evolve_ptr);
    	gsl_odeiv2_step_free(step_ptr);
    	gsl_rng_free(r);

    	MPI_Finalize();
    	
	return EXIT_SUCCESS;
}
