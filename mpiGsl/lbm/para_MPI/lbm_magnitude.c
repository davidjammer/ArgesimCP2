#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <functions.h>
#include <lbm_parameter.h>
#include <lbm_magnitude.h>

double **calc_magnitude(int Nx, int Ny, int nx, int ny, double **ux, double **uy, int dim_size[2], MPI_Comm lbm_comm)
{
	double **l_u, **u;
	int x, y, i;
	int src_rank, dest_rank, my_rank;

	int coord[2] = {0,0}; //{y,x}

	MPI_Comm_rank(lbm_comm, &my_rank);
	MPI_Cart_rank(lbm_comm, coord, &dest_rank);


	l_u = create_matrix2d(nx, ny);
	
	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			l_u[y][x] = sqrt(ux[y][x] * ux[y][x] + uy[y][x] * uy[y][x]) / u_0;
	
	u = NULL;
	
	if(dest_rank == my_rank)
	{
		u = create_matrix2d(Nx, Ny);
		
		for(coord[0]=0; coord[0]<dim_size[0]; coord[0]++)
		{
			for(coord[1]=0; coord[1]<dim_size[1]; coord[1]++)
			{
				MPI_Cart_rank(lbm_comm, coord, &src_rank);
			
				if(src_rank == my_rank)
				{
					for(i=0; i<ny; i++)
						memcpy(&u[coord[0] * ny + i][coord[1] * nx], &l_u[i][0], sizeof(double) * nx);
				}
				else
				{
					MPI_Recv(&l_u[0][0], nx * ny, MPI_DOUBLE, src_rank, 0, lbm_comm, MPI_STATUS_IGNORE);
					for(i=0; i<ny; i++)
						memcpy(&u[coord[0] * ny + i][coord[1] * nx], &l_u[i][0], sizeof(double) * nx);
				}
			}
		}
	}
	else
		MPI_Send(&l_u[0][0], nx * ny, MPI_DOUBLE, dest_rank, 0, lbm_comm);	

	free_matrix2d(l_u);

	return u;
}
