//gcc -O3 main.c `pkg-config --libs gsl`

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <string.h>

#define N_DEF 500
#define L 0.5
#define H 0.05
#define dt_DEF 0.001
#define STEPS_DEF 10000
#define NU 0.06

struct dgl
{
	double nu;
	double k;
	double h;
	int n;
};


//DGL for the 
int f(double t, const double y[], double dydt[], void *params)
{
	struct dgl *dgl;
	double nuk;
	int i;
	
	dgl = (struct dgl*) params;
	
	
	nuk = (dgl->nu * dgl->nu) / (dgl->k * dgl->k);
	
	for(i=1; i<(dgl->n); i++)
	{
		dydt[i]     = y[i + dgl->n];
		dydt[i + dgl->n] = nuk * y[i-1] - 2 * nuk * y[i]  + nuk * y[i+1]; 
	}
	
	return GSL_SUCCESS;
}

void dgl_y_init(double *y, struct dgl *dgl)
{
	int i;

	for(i=0; i<dgl->n / 2; i++)
		y[i] = 2.0 * dgl->h / (double)dgl->n * (double)i;
	for(i=dgl->n/2; i<dgl->n+1; i++)
		y[i] = 2.0 * dgl->h * (1.0 - (double)i / (double)dgl->n);

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
	int dim;
	int i;
	int steps;
	int N;
	double dt;
	int STEPS;
	struct dgl dgl;
	double *y;
	double *uxt5, *uxt8, *ut3L4, *utL2;
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
		printf("running with default parameters:\n");
		printf("N = %d\ndt = %lf\nSTEPS = %d\n", N, dt, STEPS);
	}
	
	fd=fopen("u.dat","w");
	
	
	dgl.nu = NU;
	dgl.k = L/N;
	dgl.n = N;
	dgl.h = H;
	
	//vector of the left side
	dim = 2 * (dgl.n + 1);
	y = calloc(sizeof(double), dim);
	
	//exitation at 5 and 8 seconds
	uxt5 = malloc(sizeof(double) * dim/2);
	uxt8 = malloc(sizeof(double) * dim/2);
	
	//time vector of position 3/4L and 1/2L
	ut3L4 = malloc(sizeof(double) * (STEPS+1));
	utL2 = malloc(sizeof(double) * (STEPS+1));
	
	dgl_y_init(y, &dgl);
	
	type_ptr = gsl_odeiv2_step_rk4; //RK4 Solver
	step_ptr = gsl_odeiv2_step_alloc(type_ptr, dim); //stepping function
	evolve_ptr = gsl_odeiv2_evolve_alloc(dim); //evolution function
	

	
	
	//DGL-System (jacobian is not used)
	my_system.function = f;	/* the right-hand-side functions dy[i]/dt */
    	my_system.jacobian = NULL;//jacobi;	/* the Jacobian df[i]/dy[j] */
    	my_system.dimension = dim;	/* number of diffeq's */
    	my_system.params = &dgl;	/* parameters to pass to rhs and jacobian */
	
	
	t = 0.0;
	fprintf(fd, "%E\t", 0.0);
	for(i=0; i<(dgl.n + 1); i++)
		fprintf(fd, "%E\t", (double)i * dgl.k);
	fprintf(fd,"\n");
	save(fd, y, t, dim/2);
	
	//sim loop
	for(steps=0; steps<STEPS; steps++)
	{
		gsl_odeiv2_evolve_apply_fixed_step (evolve_ptr, NULL, step_ptr,
		                            &my_system, &t, dt, y);
		
		//save the results
		save(fd, y, t, dim/2);
	}
	
	
	
	
	//final steps
	gsl_odeiv2_evolve_free (evolve_ptr);
    	gsl_odeiv2_step_free (step_ptr);
    	fclose(fd);
    	
	return EXIT_SUCCESS;
}
