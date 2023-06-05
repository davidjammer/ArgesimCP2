//gcc -O3 main.c `pkg-config --libs gsl`

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <string.h>
#include <mpi.h>

#define N_DEF 500
#define L 0.5
#define H 0.05
#define dt_DEF 0.001
#define TEND 10.0
#define STEPS_DEF 10000
#define NU 0.06


//#define DEBUG

#define RIGHT 0
#define LEFT 1

struct dgl
{
	double nu;
	double k;
	double h;
	int n;
	int lsize;
	double y_left, y_right;
};

int max(int x, int y)
{
	if(x > y)
		return x;
	else
		return y;
}

int min(int x, int y)
{
	if(x < y)
		return x;
	else
		return y;
}

//DGL for the 

int f(double t, const double y[], double dydt[], void *params)
{
	struct dgl *dgl;
	double nuk;
	
	int i;
	int world_rank, world_size;

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	dgl = (struct dgl*) params;
	
	
	nuk = (dgl->nu * dgl->nu) / (dgl->k * dgl->k);
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(world_size > 1)
	{
		//communicate the edges
		if(world_rank == 0)
		{
			MPI_Send(&y[dgl->lsize - 1], 1, MPI_DOUBLE, world_rank + 1, RIGHT, MPI_COMM_WORLD);
			MPI_Recv(&dgl->y_right, 1, MPI_DOUBLE, world_rank + 1, LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else if(world_rank == (world_size-1))
		{
			MPI_Send(&y[0], 1, MPI_DOUBLE, world_rank - 1, LEFT, MPI_COMM_WORLD);
			MPI_Recv(&dgl->y_left, 1, MPI_DOUBLE, world_rank - 1, RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Send(&y[dgl->lsize - 1], 1, MPI_DOUBLE, world_rank + 1, RIGHT, MPI_COMM_WORLD);
			MPI_Send(&y[0], 1, MPI_DOUBLE, world_rank - 1, LEFT, MPI_COMM_WORLD);
			MPI_Recv(&dgl->y_right, 1, MPI_DOUBLE, world_rank + 1, LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&dgl->y_left, 1, MPI_DOUBLE, world_rank - 1, RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		
		//calculate the right side of the DGL
		if(world_rank == 0)
		{		
			for(i=1; i<(dgl->lsize - 1); i++)
			{
				dydt[i]     = y[i + dgl->lsize];
				dydt[i + dgl->lsize] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * y[i+1]; 
			}
			dydt[i]     = y[i + dgl->lsize];
			dydt[i + dgl->lsize] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * dgl->y_right;
		}
		else if(world_rank == (world_size - 1))
		{
			i=0;
			dydt[i]     = y[i + dgl->lsize];
			dydt[i + dgl->lsize] = nuk * dgl->y_left - 2 * nuk * y[i]  + nuk * y[i+1];	
			for(i=1; i<(dgl->lsize - 1); i++)
			{
				dydt[i]     = y[i + dgl->lsize];
				dydt[i + dgl->lsize] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * y[i+1]; 
			}		
		}
		else
		{
			i=0;
			dydt[i]     = y[i + dgl->lsize];
			dydt[i + dgl->lsize] = nuk * dgl->y_left - 2 * nuk * y[i]  + nuk * y[i+1];	
			for(i=1; i<(dgl->lsize - 1); i++)
			{
				dydt[i]     = y[i + dgl->lsize];
				dydt[i + dgl->lsize] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * y[i+1]; 
			}
			dydt[i]     = y[i + dgl->lsize];
			dydt[i + dgl->lsize] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * dgl->y_right;
		}
	}
	else
	{
		for(i=1; i<(dgl->n); i++)
		{
			dydt[i]     = y[i + dgl->n];
			dydt[i + dgl->n] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * y[i+1]; 
		}
	}
	
	return GSL_SUCCESS;
}



void dgl_y_init(double *y, struct dgl *dgl)
{
	int i, imin, imax;
	int world_rank, world_size;

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	
	imin = (dgl->n + 1)  / world_size * world_rank;
	for(i=0; i<world_rank; i++)
	{
		if( ((dgl->n + 1) % world_size) > 0 && (world_size - ((dgl->n + 1) % world_size)) <= i )
			imin++;
	}
	
	
	imax = imin + dgl->lsize;

#ifdef DEBUG
	printf("%d imin: %d imax: %d\n", world_rank, imin, imax);
#endif

	for(i=max(0, imin); i<min(dgl->n / 2, imax); i++)
		y[i-imin] = 2.0 * dgl->h / (double)dgl->n * (double)i;
	for(i=max(dgl->n/2, imin); i<min(dgl->n+1, imax); i++)
		y[i-imin] = 2.0 * dgl->h * (1.0 - (double)i / (double)dgl->n);

}


//save function; save t and y in file
void save(FILE *fd, double *y, double t, int n)
{
	
	int i;
	
	fprintf(fd,"%E \t", t);
	for(i=0; i<n; i++)
	{
		fprintf(fd,"%E ",y[i]);
	}
	fprintf(fd,"\n");

}

int main(int argc, char *argv[])
{
	int world_rank, world_size;
	int dim;
	int i, n, len, imin, imax;
	int steps;
	int N;
	double dt;
	int STEPS;
	struct dgl dgl;
	double *gy, *y, *y0;
	double *ux;
	double *guxt5, *guxt8, *gut3L4, *gutL2;
	int pos_3L4, pos_L2;
	const gsl_odeiv2_step_type *type_ptr;
	gsl_odeiv2_step *step_ptr;
	gsl_odeiv2_evolve *evolve_ptr;
	gsl_odeiv2_system my_system;
	FILE *fd;
	double t;
	
	if(argc == 4)
	{
		sscanf(argv[1], "%d", &N);
		sscanf(argv[2], "%lf", &dt);
		sscanf(argv[3], "%d", &STEPS); 
	}
	else
	{
		N = N_DEF;
		dt = dt_DEF;
		STEPS = STEPS_DEF;
		
	}
	
	//Initialisierung
	MPI_Init(NULL, NULL);
    	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	if(world_rank == 0)
	{
		printf("running with default parameters:\n");
		printf("N = %d\ndt = %lf\nSTEPS = %d\n", N, dt, STEPS);
		fd=fopen("u.dat","w");
	}
	
	
	dgl.nu = NU;
	dgl.k = L / N;
	dgl.n = N;
	dgl.h = H;

	dgl.lsize = (dgl.n + 1)  / world_size;

	if( ((dgl.n + 1) % world_size) > 0 && (world_size - ((dgl.n + 1) % world_size)) <= world_rank )
		dgl.lsize++;
#ifdef DEBUG		
	printf("%3d: n=%d\n", world_rank, dgl.lsize);
#endif
	//vector of the left side
	dim = dgl.lsize * 2;
	y = calloc(sizeof(double), dim);
	
	//exitation vector
	ux = malloc(sizeof(double) * (dgl.n + 1));
	

	
	//init the vector of the left side
	dgl_y_init(y, &dgl);


		
	
	type_ptr = gsl_odeiv2_step_rk4; //RK4 Solver
	step_ptr = gsl_odeiv2_step_alloc(type_ptr, dim); //stepping function
	evolve_ptr = gsl_odeiv2_evolve_alloc(dim); //evolution function
	

	
	
	//DGL-System (jacobian is not used)
	my_system.function = f;	/* the right-hand-side functions dy[i]/dt */
    	my_system.jacobian = NULL;//jacobi;	/* the Jacobian df[i]/dy[j] */
    	my_system.dimension = dim;	/* number of diffeq's */
    	my_system.params = &dgl;	/* parameters to pass to rhs and jacobian */
	
	
	//init time vectors
	t = 0.0;
	MPI_Barrier(MPI_COMM_WORLD);
	if(world_rank == 0)
	{
		fprintf(fd, "%E\t", 0.0);
		for(i=0; i<(dgl.n + 1); i++)
			fprintf(fd, "%E\t", (double)i * dgl.k);
		fprintf(fd,"\n");
	
		memcpy(ux, y, sizeof(double) * dgl.lsize);
		len = dgl.lsize;
		for(i=1; i<world_size; i++)
		{
			MPI_Recv(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&ux[len], n, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			len+=n;
		}
		save(fd, ux, t, (dgl.n + 1));
	}
	else
	{
		MPI_Send(&dgl.lsize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(y, dgl.lsize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}
	
	
	//sim loop
	MPI_Barrier(MPI_COMM_WORLD);
	for(steps=0; steps<STEPS; steps++)
	{
	
		gsl_odeiv2_evolve_apply_fixed_step (evolve_ptr, NULL, step_ptr,
		                            &my_system, &t, dt, y);
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		if(world_rank == 0)
		{
			memcpy(ux, y, sizeof(double) * dgl.lsize);
			len = dgl.lsize;
			for(i=1; i<world_size; i++)
			{
				MPI_Recv(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&ux[len], n, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				len+=n;
			}
			save(fd, ux, t, (dgl.n + 1));
		}
		else
		{
			MPI_Send(&dgl.lsize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(y, dgl.lsize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		}
	}
	
	//final steps
	gsl_odeiv2_evolve_free (evolve_ptr);
    	gsl_odeiv2_step_free (step_ptr);
    	
    	if(world_rank == 0)
    		fclose(fd);
    	
    	MPI_Finalize();
    	
    	
	return EXIT_SUCCESS;
}


