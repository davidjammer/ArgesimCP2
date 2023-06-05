#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <functions.h>
#include <lbm_parameter.h>
#include <lbm_collision_step.h>


//calculate rho -> sum f_i [Equation 2]
void calc_rho(double **rho, double ***f, int nx, int ny)
{
	int x, y, z;

	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			rho[y][x] = 0.0;

	for(z=0; z<9; z++)
		for(y=0; y<ny; y++)
			for(x=0; x<nx; x++)
				rho[y][x] += f[z][y][x];

}

//calculate velocity in x-direction -> 1/rho * sum_i (f_i * e_i) [Equation 3]
void calc_ux(double **ux, double ***f, double **rho, int nx, int ny)
{
	int x, y;
	
	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			ux[y][x] = (f[E][y][x] - f[W][y][x] + f[NE][y][x] + f[SE][y][x] - f[SW][y][x] - f[NW][y][x]) / rho[y][x];

}

//calculate velocity in y-direction -> 1/rho * sum_i (f_i * e_i) [Equation 3]
void calc_uy(double **uy, double ***f, double **rho, int nx, int ny)
{
	int x, y;
	
	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			uy[y][x] = (f[N][y][x] - f[S][y][x] + f[NE][y][x] + f[NW][y][x] - f[SE][y][x] - f[SW][y][x]) / rho[y][x];

}

//set the velocity from the driving cells in x-direction to u_0
void set_ux_for_driving_cells(double **ux, index2d *driving_ids, int n_driving_ids)
{
	int i;
	
	for(i=0; i<n_driving_ids; i++)
		ux[driving_ids[i].y][driving_ids[i].x] = u_0;
	
}

//set the velocity from the driving cells in y-direction to 0
void set_uy_for_driving_cells(double **uy, index2d *driving_ids, int n_driving_ids)
{
	int i;
	
	for(i=0; i<n_driving_ids; i++)
		uy[driving_ids[i].y][driving_ids[i].x] = 0.0;

}

//calculate the helper matrix usqr: ux^2 + uy^2
void calc_usqr(double **usqr, double **ux, double **uy, int nx, int ny)
{
	int x, y;

	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			usqr[y][x] = ux[y][x] * ux[y][x] + uy[y][x] * uy[y][x];

}

//calulate the equilibrium distribution [Equation 4]
void clac_feq(double ***feq, double **rho, double **ux, double **uy, double **usqr, int nx, int ny)
{
	int x, y;

	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			feq[C][y][x] = 4.0 / 9.0 * rho[y][x] * (1.0 - 1.5 * usqr[y][x]);

	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			feq[E][y][x] = 1.0 / 9.0 * rho[y][x] * (1.0 + 3.0 * ux[y][x] + 4.5 * ux[y][x] * ux[y][x] - 1.5 * usqr[y][x]);
	
	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			feq[S][y][x] = 1.0 / 9.0 * rho[y][x] * (1.0 - 3.0 * uy[y][x] + 4.5 * uy[y][x] * uy[y][x] - 1.5 * usqr[y][x]);

	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			feq[W][y][x] = 1.0 / 9.0 * rho[y][x] * (1.0 - 3.0 * ux[y][x] + 4.5 * ux[y][x] * ux[y][x] - 1.5 * usqr[y][x]);

	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			feq[N][y][x] = 1.0 / 9.0 * rho[y][x] * (1.0 + 3.0 * uy[y][x] + 4.5 * uy[y][x] * uy[y][x] - 1.5 * usqr[y][x]);

	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			feq[NE][y][x] = 1.0 / 36.0 * rho[y][x] * (1.0 + 3.0 * (ux[y][x] + uy[y][x]) + 4.5 * (ux[y][x] + uy[y][x]) * (ux[y][x] + uy[y][x]) - 1.5 * usqr[y][x]);
		
	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			feq[SE][y][x] = 1.0 / 36.0 * rho[y][x] * (1.0 + 3.0 * (ux[y][x] - uy[y][x]) + 4.5 * (ux[y][x] - uy[y][x]) * (ux[y][x] - uy[y][x]) - 1.5 * usqr[y][x]);

	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			feq[SW][y][x] = 1.0 / 36.0 * rho[y][x] * (1.0 + 3.0 * (-ux[y][x] - uy[y][x]) + 4.5 * (-ux[y][x] - uy[y][x]) * (-ux[y][x] - uy[y][x]) - 1.5 * usqr[y][x]);

	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			feq[NW][y][x] = 1.0 / 36.0 * rho[y][x] * (1.0 + 3.0 * (-ux[y][x] + uy[y][x]) + 4.5 * (-ux[y][x] + uy[y][x]) * (-ux[y][x] + uy[y][x]) - 1.5 * usqr[y][x]);
}

//calculate the distribution value according to cell types
//wall cells -> bounce back
//driving cells -> equilibrium values
//fluid cells -> [Equation 1]
void f_cell_calc(double ***f, double ***feq, double **geometry, double tau, int nx, int ny)
{
	double temp;
	int x, y;
	
	for(y=0; y<ny; y++)
	{
		for(x=0; x<nx; x++)
		{
			if(geometry[y][x] == WALL_CELL) //Wall Cell
			{
				temp = f[E][y][x];
				f[E][y][x] = f[W][y][x];
				f[W][y][x] = temp;
				
				temp = f[S][y][x];
				f[S][y][x] = f[N][y][x];
				f[N][y][x] = temp;
				
				temp = f[NE][y][x];
				f[NE][y][x] = f[SW][y][x];
				f[SW][y][x] = temp;
				
				temp = f[NW][y][x];
				f[NW][y][x] = f[SE][y][x];
				f[SE][y][x] = temp;
			}	
			else if(geometry[y][x] == DRIVING_CELL) //Driving Cell
			{
				f[C][y][x] = feq[C][y][x];
				f[N][y][x] = feq[N][y][x];
				f[E][y][x] = feq[E][y][x];
				f[S][y][x] = feq[S][y][x];
				f[W][y][x] = feq[W][y][x];
				f[NE][y][x] = feq[NE][y][x];
				f[SE][y][x] = feq[SE][y][x];
				f[SW][y][x] = feq[SW][y][x];
				f[NW][y][x] = feq[NW][y][x];
			}
			else //Fluid Cell
			{
				f[C][y][x] = f[C][y][x] * (1.0 - 1.0 / tau) + feq[C][y][x] / tau;
				
				f[N][y][x] = f[N][y][x] * (1.0 - 1.0 / tau) + feq[N][y][x] / tau;
				f[E][y][x] = f[E][y][x] * (1.0 - 1.0 / tau) + feq[E][y][x] / tau;
				f[S][y][x] = f[S][y][x] * (1.0 - 1.0 / tau) + feq[S][y][x] / tau;
				f[W][y][x] = f[W][y][x] * (1.0 - 1.0 / tau) + feq[W][y][x] / tau;
				
				f[NE][y][x] = f[NE][y][x] * (1.0 - 1.0 / tau) + feq[NE][y][x] / tau;
				f[SE][y][x] = f[SE][y][x] * (1.0 - 1.0 / tau) + feq[SE][y][x] / tau;
				f[SW][y][x] = f[SW][y][x] * (1.0 - 1.0 / tau) + feq[SW][y][x] / tau;
				f[NW][y][x] = f[NW][y][x] * (1.0 - 1.0 / tau) + feq[NW][y][x] / tau;
			}
		}
	}
}
