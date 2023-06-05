#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <functions.h>
#include <lbm_parameter.h>
#include <lbm_propagation_step.h>

//shift the matrix m one step to the north-east
void shift_NE(double **m, int nx, int ny)
{
	int y;
	
	for(y = 1; y<ny; y++)
	{
		memcpy(&m[y-1][1], &m[y][0], sizeof(double) * (nx - 1));
	}
}

//shift the matrix m one step to the south-east
void shift_SE(double **m, int nx, int ny)
{
	int y;
	
	for(y = ny-2; y>=0; y--)
	{
		memcpy(&m[y+1][1], &m[y][0], sizeof(double) * (nx - 1));
	}
}

//shift the matrix m one step to the south-west
void shift_SW(double **m, int nx, int ny)
{
	int y;
	
	for(y = ny-2; y>=0; y--)
	{
		memcpy(&m[y+1][0], &m[y][1], sizeof(double) * (nx - 1));
	}
}

//shift the matrix m one step to the north-west
void shift_NW(double **m, int nx, int ny)
{
	int y;
	
	for(y = 0; y<(ny-1); y++)
	{
		memcpy(&m[y][0], &m[y+1][1], sizeof(double) * (nx - 1));
	}
}

//shift the matrix m one step to the east
void shift_E(double **m, int nx, int ny)
{
	int y;
	double *temp;
	
	temp = malloc(sizeof(double) * (nx-1));
	
	for(y=0; y<ny; y++)
	{
		memcpy(temp, &m[y][0], sizeof(double) * (nx-1));
		memcpy(&m[y][1], temp, sizeof(double) * (nx-1));
	}
	
	free(temp);
}

//shift the matrix m one step to the west
void shift_W(double **m, int nx, int ny)
{
	int y;
	double *temp;
	
	temp = malloc(sizeof(double) * (nx-1));
	
	for(y=0; y<ny; y++)
	{
		memcpy(temp, &m[y][1], sizeof(double) * (nx-1));
		memcpy(&m[y][0], temp, sizeof(double) * (nx-1));
	}
	
	free(temp);
}

//shift the matrix m one step to the south
void shift_S(double **m, int nx, int ny)
{
	int y;
	
	for(y=(ny-2); y>=0; y--)
	{
		memcpy(&m[y+1][0], &m[y][0], sizeof(double) * (nx));
	}
}

//shift the matrix m one step to the north 
void shift_N(double **m, int nx, int ny)
{
	int y;
	
	for(y=0; y<(ny-1); y++)
	{
		memcpy(&m[y][0], &m[y+1][0], sizeof(double) * (nx));
	}
}
/*do the propagation step
*particle moving:
* E
* S
* W
* N
* NE
* SE
* SW
* NW
*/
void propagation_step(double ***f, int nx, int ny)
{
	//East	
	shift_E(f[E], nx, ny);
	//South
	shift_S(f[S], nx, ny);
	//West
	shift_W(f[W], nx, ny);
	//North
	shift_N(f[N], nx, ny);
	//North-East
	shift_NE(f[NE], nx, ny);
	//South-East
	shift_SE(f[SE], nx, ny);
	//South-West
	shift_SW(f[SW], nx, ny);
	//North-West
	shift_NW(f[NW], nx, ny);
}
