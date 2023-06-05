#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <functions.h>
#include <lbm_init.h>
#include <lbm_parameter.h>

/*create the geometry matrix (g) for the benchmark
*
* 2 -> driving cell 
* 1 -> wall cell
* 0 -> fluid cell
*
* 2 2 2 2 2 2 
* 1 0 0 0 0 1
* 1 0 0 0 0 1
* 1 0 0 0 0 1
* 1 0 0 0 0 1
* 1 1 1 1 1 1
*/
void create_geometry(double **g, int nx, int ny)
{
	int x, y;

	//wall cells
	for(x=0; x<nx; x++)
		g[ny-1][x] = WALL_CELL;
	for(y=1; y<(ny-1); y++)
		g[y][0] = WALL_CELL;
	for(y=1; y<(ny-1); y++)
		g[y][nx-1] = WALL_CELL;

	//fluid cells
	for(y=1; y<(ny-1); y++)
		for(x=1; x<(nx-1); x++)
			g[y][x] = FLUID_CELL;
	
	//driving cells
	for(x=0; x<nx; x++)
		g[0][x] = DRIVING_CELL;


}

/*initializes the matrix f
* non moving particles [C]: 4/9
* moving particles [E S W N]: 1/9
* moving particles diagonal [NE SE SW NE]: 1/36
* Paper Rene: 
* non moving [C] -> a:1/9 
* moving particle [E S W N] -> a:4/9 
* moving particle [NE SE SW NW] -> a:1/36
*/
void init_dist_func_values(double ***f, int nx, int ny)
{
	int x, y;
	
	set_value_on_matrix(f[C], nx, ny, RHO_0 * 4.0 / 9.0);
	
	set_value_on_matrix(f[E], nx, ny, RHO_0 / 9.0);
	set_value_on_matrix(f[S], nx, ny, RHO_0 / 9.0);
	set_value_on_matrix(f[W], nx, ny, RHO_0 / 9.0);
	set_value_on_matrix(f[N], nx, ny, RHO_0 / 9.0);
	
	set_value_on_matrix(f[NE], nx, ny, RHO_0 / 36.0);
	set_value_on_matrix(f[SE], nx, ny, RHO_0 / 36.0);
	set_value_on_matrix(f[SW], nx, ny, RHO_0 / 36.0);
	set_value_on_matrix(f[NW], nx, ny, RHO_0 / 36.0);
}
