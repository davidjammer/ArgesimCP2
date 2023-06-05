#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

#include <functions.h>
#include <lbm_propagation_step.h>

#define Nx 8
#define Ny 8

int main(int argc, char *argv[])
{
	MPI_Comm lbm_comm;
	int world_size, world_rank, my_rank, src_rank, dest_rank;
	double **g_geometry, **l_geometry,**u;
	int x,y,i;
	int nx, ny;
	
	int dim_size[2] = {0,0};
	int periods[2] = {0,0};
	int coord[2] = {0,0}; //{y,x}

	//Initialization of MPI
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);



	MPI_Dims_create(world_size, 2, dim_size);
	
	MPI_Cart_create(MPI_COMM_WORLD,2,dim_size,periods,0,&lbm_comm);
	
	ny = Ny / dim_size[0];
	nx = Nx / dim_size[1];

	MPI_Comm_rank(lbm_comm, &my_rank);
	MPI_Cart_rank(lbm_comm, coord, &src_rank);

	l_geometry = create_matrix2d(nx, ny);
	if(src_rank == my_rank)
	{
		g_geometry = create_matrix2d(Nx, Ny);
		
		for(y=0, i=1; y<Ny; y++)
		 	for(x=0; x<Nx; x++)
				g_geometry[y][x] = i++;
		
		for(coord[0]=dim_size[0]-1; coord[0]>=0; coord[0]--)
		{
			for(coord[1]=dim_size[1]-1; coord[1]>=0; coord[1]--)
			{
				//copy matrix
				for(i=0; i<ny; i++)
					memcpy(l_geometry[i], &g_geometry[coord[0] * ny + i][coord[1] * nx], sizeof(double) * nx);
				if(!(coord[0]==0 && coord[1]==0))
				{ 
					MPI_Cart_rank(lbm_comm, coord, &dest_rank);
					MPI_Send(&l_geometry[0][0], nx * ny, MPI_DOUBLE, dest_rank, 0, lbm_comm);
				}
			}
		}
	}
	else
		MPI_Recv(&l_geometry[0][0], nx * ny, MPI_DOUBLE, 0, src_rank, lbm_comm, MPI_STATUS_IGNORE);
	
	//shift_N(l_geometry, nx, ny, lbm_comm);
	//shift_E(l_geometry, nx, ny, lbm_comm);
	//shift_S(l_geometry, nx, ny, lbm_comm);
	//shift_W(l_geometry, nx, ny, lbm_comm);
	
	//shift_NE(l_geometry, nx, ny, lbm_comm);
	//shift_SE(l_geometry, nx, ny, lbm_comm);
	//shift_SW(l_geometry, nx, ny, lbm_comm);
	shift_NW(l_geometry, nx, ny, lbm_comm);
	
	
	coord[0] = 0;
	coord[1] = 0;

	MPI_Comm_rank(lbm_comm, &my_rank);
	MPI_Cart_rank(lbm_comm, coord, &dest_rank);
	
	u = NULL;
	
	printf("dest_rank %d my_rank %d\n",dest_rank,my_rank);
	
	if(dest_rank == my_rank)
	{
		u = create_matrix2d(Nx, Ny);
		printf("create u\n");
		for(coord[0]=0; coord[0]<dim_size[0]; coord[0]++)
		{
			for(coord[1]=0; coord[1]<dim_size[1]; coord[1]++)
			{
				MPI_Cart_rank(lbm_comm, coord, &src_rank);
			
				if(src_rank == my_rank)
				{
					for(i=0; i<ny; i++)
						memcpy(&u[coord[0] * ny + i][coord[1] * nx], &l_geometry[i][0], sizeof(double) * nx);
				}
				else
				{
					MPI_Recv(&l_geometry[0][0], nx * ny, MPI_DOUBLE, src_rank, 0, lbm_comm, MPI_STATUS_IGNORE);
					for(i=0; i<ny; i++)
						memcpy(&u[coord[0] * ny + i][coord[1] * nx], &l_geometry[i][0], sizeof(double) * nx);
				}
			}
		}
	}
	else
		MPI_Send(&l_geometry[0][0], nx * ny, MPI_DOUBLE, dest_rank, 0, lbm_comm);
	
	if(u != NULL)
		save("ushift.dat", u, Nx, Ny);
		
		
	MPI_Finalize();
	
	return EXIT_SUCCESS;
	
}
