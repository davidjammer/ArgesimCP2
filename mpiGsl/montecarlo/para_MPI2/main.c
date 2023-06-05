//mpicc -O3 main.c `pkg-config --libs --cflags gsl`

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <time.h>
#include <mpi.h>

#define TIME 60.0   //max. simulation time
#define STEPS 200 //number of compute steps
#define H 0.01   //stepsize

//parameter for the damped dynamic mass - spring system
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

//jacobian matrix for the damped dynamic mass - spring system
int jacobi(double t, const double y[], double * dfdy, double dfdt[], void * params)
{
	double k, d, m;
	gsl_matrix *m_ptr;
	gsl_matrix_view dfdy_mat;
	
	k = ((double*)params)[0];
	d = ((double*)params)[1];
	m = ((double*)params)[2];
	
	dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);

    	m_ptr = &dfdy_mat.matrix;	/* m_ptr points to the matrix */

    	/* fill the Jacobian matrix as shown */
    	gsl_matrix_set (m_ptr, 0, 0, 0.0);	/* df[0]/dy[0] = 0 */
    	gsl_matrix_set (m_ptr, 0, 1, 1.0);	/* df[0]/dy[1] = 1 */
    	gsl_matrix_set (m_ptr, 1, 0, - k / m); /* df[1]/dy[0] */
    	gsl_matrix_set (m_ptr, 1, 1, - d / m);     /* df[1]/dy[1] */

    	/* set explicit t dependence of f[i] */
    	dfdt[0] = 0.0;
 	dfdt[1] = 0.0;


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

double calc_t(struct timespec t_start, struct timespec t)
{
	double t1, t2;
	
	t1 = (double)t.tv_sec;
	t2 = (double)t_start.tv_sec;
	
	return t1 - t2;
}

int main()
{
	int world_size;
	int world_rank;
	int dim = 2;
	int i;
	int steps;
	double *tvec;
	double *y0vec, *y0vec_recv;
	const gsl_odeiv2_step_type *type_ptr;
	gsl_odeiv2_step *step_ptr;
	gsl_odeiv2_evolve *evolve_ptr;
	gsl_odeiv2_system my_system;
	double y[2];
	double t;
	double params[3];
	struct timespec t_start, t_now;
	unsigned int num_sim, sum_sim;
	 
	
	
	//Initialisierung
	MPI_Init(NULL, NULL);
    	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	//vector for the time and the mean motion
	tvec = calloc(sizeof(double), STEPS);
	y0vec = calloc(sizeof(double), STEPS);

	
	type_ptr = gsl_odeiv2_step_rk4; //RK4 Solver
	step_ptr = gsl_odeiv2_step_alloc(type_ptr, dim); //stepping function
	evolve_ptr = gsl_odeiv2_evolve_alloc(dim); //evolution function
	

	//parameter vector for the DGL
	params[0] = K;
	params[2] = M; 
	
	//DGL-System (jacobian is not used)
	my_system.function = f;	/* the right-hand-side functions dy[i]/dt */
    	my_system.jacobian = NULL;//jacobi;	/* the Jacobian df[i]/dy[j] */
    	my_system.dimension = dim;	/* number of diffeq's */
    	my_system.params = params;	/* parameters to pass to rhs and jacobian */
	
	

	clock_gettime(CLOCK_MONOTONIC, &t_start);
	num_sim = 0;
	do
	{
		//set the state vector motion=0 and velocity=0.1 
		y[0] = 0.0;
		y[1] = 0.1;
		//set time
		t = 0.0;
	
		//set damping factor
		params[1] = (D_MAX - D_MIN) * ((double)rand() / (double)RAND_MAX) + D_MIN;
		
		//first values
		tvec[0] = t;
		y0vec[0] += y[0];
		
		//sim. loop
		for(steps=1; steps < STEPS; steps++)
		{
			//ODE-Solver
			gsl_odeiv2_evolve_apply_fixed_step (evolve_ptr, NULL, step_ptr,
		                            &my_system, &t, H, y);

			//set result for time t
		    	tvec[steps] = t;
			y0vec[steps] += y[0];
		}
		num_sim++;
		clock_gettime(CLOCK_MONOTONIC, &t_now);
	}while(calc_t(t_start, t_now) <= TIME);
	MPI_Barrier(MPI_COMM_WORLD);

	printf("%d -> %d simulations\n", world_rank, num_sim);
	if(world_rank == 0)	
		MPI_Reduce(&num_sim, &sum_sim, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
	else
		MPI_Reduce(&num_sim, NULL, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if(world_rank == 0)
		printf("Sum.: %d simulations\n", sum_sim);

	if(world_rank == 0)
	{
		y0vec_recv = malloc(sizeof(double) * STEPS);
		MPI_Reduce(y0vec, y0vec_recv, STEPS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		//compute the mean motion
		for(steps=0; steps < STEPS; steps++)
			y0vec_recv[steps] /= (double)sum_sim;

		//save results
		save("daten.dat", tvec, y0vec_recv, STEPS);
	}
	else
		MPI_Reduce(y0vec, NULL, STEPS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	
		
	//final steps
	gsl_odeiv2_evolve_free (evolve_ptr);
    	gsl_odeiv2_step_free (step_ptr);
    	
    	MPI_Finalize();
    	
	return EXIT_SUCCESS;
}
