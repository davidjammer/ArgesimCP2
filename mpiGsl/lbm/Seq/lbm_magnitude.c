#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <functions.h>
#include <lbm_parameter.h>
#include <lbm_magnitude.h>

double **calc_magnitude(int Nx, int Ny, double **ux, double **uy)
{
	double **u;
	int x, y;

	u = create_matrix2d(Nx, Ny);
	
	for(y=0; y<Nx; y++)
		for(x=0; x<Ny; x++)
			u[y][x] = sqrt(ux[y][x] * ux[y][x] + uy[y][x] * uy[y][x]) / u_0;
			
	return u;
}
