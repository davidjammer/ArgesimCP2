#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <functions.h>
#include <lbm_parameter.h>
#include <lbm_propagation_step.h>

//shift the matrix m one step to the north-east
//send the north border, east border and the north-east point to the neighbors
//receive the south border, west border and the south-west point from the neighbors
void shift_NE(double **m, int nx, int ny, MPI_Comm lbm_comm)
{
	double *north_border, *east_border, ne_point, *data;
	int y;
	int dims[2], periods[2], coords[2];
	int my_coord[2], dest_coord[2], src_coord[2], my_rank;
	int src_S, dest_N, src_W, dest_E, src_SW, dest_NE;
	
	MPI_Cart_get(lbm_comm, 2, dims, periods, coords);
	
	north_border = malloc(sizeof(double) * (nx - 1));
	east_border = malloc(sizeof(double) * (ny - 1));
	data = malloc(sizeof(double) * (max(nx,ny) - 1));
	
	memcpy(north_border, &m[0][0], sizeof(double) * (nx - 1));
	ne_point = m[0][nx-1];

	for(y = 1; y<ny; y++)
	{
		east_border[y-1] = m[y][nx-1];
		memcpy(&m[y-1][1], &m[y][0], sizeof(double) * (nx - 1));
	}

	MPI_Cart_shift(lbm_comm, 0, -1, &src_S, &dest_N);
	MPI_Cart_shift(lbm_comm, 1, 1, &src_W, &dest_E);
	
	MPI_Comm_rank(lbm_comm, &my_rank);
	MPI_Cart_coords(lbm_comm, my_rank, 2, my_coord);
	dest_coord[0] = my_coord[0] - 1;
	dest_coord[1] = my_coord[1] + 1;
	src_coord[0] = my_coord[0] + 1;
	src_coord[1] = my_coord[1] - 1;
	
	if(dest_coord[0] >= 0 && dest_coord[0] < dims[0] && dest_coord[1] >= 0 && dest_coord[1] < dims[1])
		MPI_Cart_rank(lbm_comm, dest_coord, &dest_NE);
	else
		dest_NE = -1;
	if(src_coord[0] >= 0 && src_coord[0] < dims[0] && src_coord[1] >= 0 && src_coord[1] < dims[1])	
		MPI_Cart_rank(lbm_comm, src_coord, &src_SW);
	else
		src_SW = -1;
	
		
	if(dest_N >= 0)
	{
		MPI_Send(north_border, nx - 1, MPI_DOUBLE, dest_N, N, lbm_comm);
	}
	if(src_S >= 0)
	{

		MPI_Recv(data, nx - 1, MPI_DOUBLE, src_S, N, lbm_comm, MPI_STATUS_IGNORE);
		memcpy(&m[ny-1][1], data, sizeof(double) * (nx - 1));
	}
	
	if(dest_E >= 0)
	{
		MPI_Send(east_border, ny - 1, MPI_DOUBLE, dest_E, E, lbm_comm);
	}
	if(src_W >= 0)
	{

		MPI_Recv(data, ny - 1, MPI_DOUBLE, src_W, E, lbm_comm, MPI_STATUS_IGNORE);
		for(y=0; y<ny-1; y++)
			m[y][0] = data[y];
	}

	if(dest_NE >= 0)
	{
		MPI_Send(&ne_point, 1, MPI_DOUBLE, dest_NE, N+E, lbm_comm);
	}
	if(src_SW >= 0)
	{
		MPI_Recv(&m[ny-1][0], 1, MPI_DOUBLE, src_SW, N+E, lbm_comm, MPI_STATUS_IGNORE);
	}
			
	free(north_border);
	free(east_border);
	free(data);
	
}

//shift the matrix m one step to the south-east
//send the south border, east border and the south-east point to the neighbors
//receive the north border, west border and the north-west point from the neighbors
void shift_SE(double **m, int nx, int ny, MPI_Comm lbm_comm)
{
	double *south_border, *east_border, se_point, *data;
	int y;
	int dims[2], periods[2], coords[2];
	int my_coord[2], dest_coord[2], src_coord[2], my_rank;
	int src_N, dest_S, src_W, dest_E, src_NW, dest_SE;
	
	MPI_Cart_get(lbm_comm, 2, dims, periods, coords);
	
	south_border = malloc(sizeof(double) * (nx - 1));
	east_border = malloc(sizeof(double) * (ny - 1));
	data = malloc(sizeof(double) * (max(nx,ny) - 1));
	
	memcpy(south_border, &m[ny-1][0], sizeof(double) * (nx - 1));
	se_point = m[ny-1][nx-1];

	for(y = ny-2; y>=0; y--)
	{
		east_border[y] = m[y][nx-1];
		memcpy(&m[y+1][1], &m[y][0], sizeof(double) * (nx - 1));
	}

	MPI_Cart_shift(lbm_comm, 0, 1, &src_N, &dest_S);
	MPI_Cart_shift(lbm_comm, 1, 1, &src_W, &dest_E);
	
	MPI_Comm_rank(lbm_comm, &my_rank);
	MPI_Cart_coords(lbm_comm, my_rank, 2, my_coord);
	dest_coord[0] = my_coord[0] + 1;
	dest_coord[1] = my_coord[1] + 1;
	src_coord[0] = my_coord[0] - 1;
	src_coord[1] = my_coord[1] - 1;
	
	if(dest_coord[0] >= 0 && dest_coord[0] < dims[0] && dest_coord[1] >= 0 && dest_coord[1] < dims[1])
		MPI_Cart_rank(lbm_comm, dest_coord, &dest_SE);
	else
		dest_SE = -1;
	if(src_coord[0] >= 0 && src_coord[0] < dims[0] && src_coord[1] >= 0 && src_coord[1] < dims[1])	
		MPI_Cart_rank(lbm_comm, src_coord, &src_NW);
	else
		src_NW = -1;
	
		
	if(dest_S >= 0)
	{
		MPI_Send(south_border, nx - 1, MPI_DOUBLE, dest_S, S, lbm_comm);
	}
	if(src_N >= 0)
	{

		MPI_Recv(data, nx - 1, MPI_DOUBLE, src_N, S, lbm_comm, MPI_STATUS_IGNORE);
		memcpy(&m[0][1], data, sizeof(double) * (nx - 1));
	}
	
	if(dest_E >= 0)
	{
		MPI_Send(east_border, ny - 1, MPI_DOUBLE, dest_E, E, lbm_comm);
	}
	if(src_W >= 0)
	{

		MPI_Recv(data, ny - 1, MPI_DOUBLE, src_W, E, lbm_comm, MPI_STATUS_IGNORE);
		for(y=1; y<ny; y++)
			m[y][0] = data[y-1];
	}

	if(dest_SE >= 0)
	{
		MPI_Send(&se_point, 1, MPI_DOUBLE, dest_SE, S+E, lbm_comm);
	}
	if(src_NW >= 0)
	{
		MPI_Recv(&m[0][0], 1, MPI_DOUBLE, src_NW, S+E, lbm_comm, MPI_STATUS_IGNORE);
	}
		
	free(south_border);
	free(east_border);
	free(data);
	
}

//shift the matrix m one step to the south-west
//send the south border, west border and the south-west point to the neighbors
//receive the north border, east border and the north-east point from the neighbors
void shift_SW(double **m, int nx, int ny, MPI_Comm lbm_comm)
{
	double *south_border, *west_border, sw_point, *data;
	int y;
	int dims[2], periods[2], coords[2];
	int my_coord[2], dest_coord[2], src_coord[2], my_rank;
	int src_N, dest_S, src_E, dest_W, src_NE, dest_SW;
	
	MPI_Cart_get(lbm_comm, 2, dims, periods, coords);
	
	south_border = malloc(sizeof(double) * (nx - 1));
	west_border = malloc(sizeof(double) * (ny - 1));
	data = malloc(sizeof(double) * (max(nx,ny) - 1));
	
	memcpy(south_border, &m[ny-1][1], sizeof(double) * (nx - 1));
	sw_point = m[ny-1][0];

	for(y = ny-2; y>=0; y--)
	{
		west_border[y] = m[y][0];
		memcpy(&m[y+1][0], &m[y][1], sizeof(double) * (nx - 1));
	}

	MPI_Cart_shift(lbm_comm, 0, 1, &src_N, &dest_S);
	MPI_Cart_shift(lbm_comm, 1, -1, &src_E, &dest_W);
	
	MPI_Comm_rank(lbm_comm, &my_rank);
	MPI_Cart_coords(lbm_comm, my_rank, 2, my_coord);
	dest_coord[0] = my_coord[0] + 1;
	dest_coord[1] = my_coord[1] - 1;
	src_coord[0] = my_coord[0] - 1;
	src_coord[1] = my_coord[1] + 1;
	
	if(dest_coord[0] >= 0 && dest_coord[0] < dims[0] && dest_coord[1] >= 0 && dest_coord[1] < dims[1])
		MPI_Cart_rank(lbm_comm, dest_coord, &dest_SW);
	else
		dest_SW = -1;
	if(src_coord[0] >= 0 && src_coord[0] < dims[0] && src_coord[1] >= 0 && src_coord[1] < dims[1])	
		MPI_Cart_rank(lbm_comm, src_coord, &src_NE);
	else
		src_NE = -1;
	
		
	if(dest_S >= 0)
	{
		MPI_Send(south_border, nx - 1, MPI_DOUBLE, dest_S, S, lbm_comm);
	}
	if(src_N >= 0)
	{

		MPI_Recv(data, nx - 1, MPI_DOUBLE, src_N, S, lbm_comm, MPI_STATUS_IGNORE);
		memcpy(&m[0][0], data, sizeof(double) * (nx - 1));
	}
	
	if(dest_W >= 0)
	{
		MPI_Send(west_border, ny - 1, MPI_DOUBLE, dest_W, W, lbm_comm);
	}
	if(src_E >= 0)
	{

		MPI_Recv(data, ny - 1, MPI_DOUBLE, src_E, W, lbm_comm, MPI_STATUS_IGNORE);
		for(y=1; y<ny; y++)
			m[y][nx-1] = data[y-1];
	}

	if(dest_SW >= 0)
	{
		MPI_Send(&sw_point, 1, MPI_DOUBLE, dest_SW, S+W, lbm_comm);
	}
	if(src_NE >= 0)
	{
		MPI_Recv(&m[0][nx-1], 1, MPI_DOUBLE, src_NE, S+W, lbm_comm, MPI_STATUS_IGNORE);
	}
	
	free(south_border);
	free(west_border);
	free(data);
	
}

//shift the matrix m one step to the north-west
//send the north border, west border and the north-west point to the neighbors
//receive the south border, east border and the south-east point from the neighbors
void shift_NW(double **m, int nx, int ny, MPI_Comm lbm_comm)
{
	double *north_border, *west_border, nw_point, *data;
	int y;
	int dims[2], periods[2], coords[2];
	int my_coord[2], dest_coord[2], src_coord[2], my_rank;
	int src_S, dest_N, src_E, dest_W, src_SE, dest_NW;
	
	MPI_Cart_get(lbm_comm, 2, dims, periods, coords);
	
	north_border = malloc(sizeof(double) * (nx - 1));
	west_border = malloc(sizeof(double) * (ny - 1));
	data = malloc(sizeof(double) * (max(nx,ny) - 1));
	
	memcpy(north_border, &m[0][1], sizeof(double) * (nx - 1));
	nw_point = m[0][0];

	for(y = 0; y<(ny-1); y++)
	{
		west_border[y] = m[y+1][0];
		memcpy(&m[y][0], &m[y+1][1], sizeof(double) * (nx - 1));
	}

	MPI_Cart_shift(lbm_comm, 0, -1, &src_S, &dest_N);
	MPI_Cart_shift(lbm_comm, 1, -1, &src_E, &dest_W);
	
	MPI_Comm_rank(lbm_comm, &my_rank);
	MPI_Cart_coords(lbm_comm, my_rank, 2, my_coord);
	dest_coord[0] = my_coord[0] - 1;
	dest_coord[1] = my_coord[1] - 1;
	src_coord[0] = my_coord[0] + 1;
	src_coord[1] = my_coord[1] + 1;
	
	if(dest_coord[0] >= 0 && dest_coord[0] < dims[0] && dest_coord[1] >= 0 && dest_coord[1] < dims[1])
		MPI_Cart_rank(lbm_comm, dest_coord, &dest_NW);
	else
		dest_NW = -1;
	if(src_coord[0] >= 0 && src_coord[0] < dims[0] && src_coord[1] >= 0 && src_coord[1] < dims[1])	
		MPI_Cart_rank(lbm_comm, src_coord, &src_SE);
	else
		src_SE = -1;
	
		
	if(dest_N >= 0)
	{
		MPI_Send(north_border, nx - 1, MPI_DOUBLE, dest_N, N, lbm_comm);
	}
	if(src_S >= 0)
	{
		MPI_Recv(data, nx - 1, MPI_DOUBLE, src_S, N, lbm_comm, MPI_STATUS_IGNORE);
		memcpy(&m[ny-1][0], data, sizeof(double) * (nx - 1));
	}
	
	if(dest_W >= 0)
	{
		MPI_Send(west_border, ny - 1, MPI_DOUBLE, dest_W, W, lbm_comm);
	}
	if(src_E >= 0)
	{
		MPI_Recv(data, ny - 1, MPI_DOUBLE, src_E, W, lbm_comm, MPI_STATUS_IGNORE);
		for(y=0; y<(ny-1); y++)
			m[y][nx-1] = data[y];
	}

	if(dest_NW >= 0)
	{
		MPI_Send(&nw_point, 1, MPI_DOUBLE, dest_NW, N+W, lbm_comm);
	}
	if(src_SE >= 0)
	{
		MPI_Recv(&m[ny-1][nx-1], 1, MPI_DOUBLE, src_SE, N+W, lbm_comm, MPI_STATUS_IGNORE);
	}
	
	free(north_border);
	free(west_border);
	free(data);
}

//shift the matrix m one step to the east
//send the east border to the neighbor
//receive the west border from the neighbor
void shift_E(double **m, int nx, int ny, MPI_Comm lbm_comm)
{
	int y;
	int src_W, dest_E;
	double *temp;
	double *data;
	
	
	MPI_Cart_shift(lbm_comm, 1, 1, &src_W, &dest_E);
	
	
	temp = malloc(sizeof(double) * (nx-1));
	data  = malloc(sizeof(double) * (ny));
	
	//shift matrix
	for(y=0; y<ny; y++)
	{
		data[y] = m[y][nx-1]; //right border
		memcpy(temp, &m[y][0], sizeof(double) * (nx-1));
		memcpy(&m[y][1], temp, sizeof(double) * (nx-1));
	}
	
	free(temp);
	
	if(dest_E >= 0)
	{
		MPI_Send(data, ny, MPI_DOUBLE, dest_E, E, lbm_comm);
	}
	if(src_W >= 0)
	{

		MPI_Recv(data, ny, MPI_DOUBLE, src_W, E, lbm_comm, MPI_STATUS_IGNORE);
		for(y=0; y<ny; y++)
			m[y][0] = data[y];
	}
	
	free(data);
}

//shift the matrix m one step to the west
//send the west border to the neighbor
//receive the east border from the neighbor
void shift_W(double **m, int nx, int ny, MPI_Comm lbm_comm)
{
	int y;
	int src_E, dest_W;
	double *temp;
	double *data;
	
	MPI_Cart_shift(lbm_comm, 1, -1, &src_E, &dest_W);
	
	temp = malloc(sizeof(double) * (nx-1));
	data  = malloc(sizeof(double) * (ny));
	
	for(y=0; y<ny; y++)
	{
		data[y] = m[y][0];
		memcpy(temp, &m[y][1], sizeof(double) * (nx-1));
		memcpy(&m[y][0], temp, sizeof(double) * (nx-1));
	}
	
	free(temp);
	
	if(dest_W >= 0)
	{	
		MPI_Send(data, ny, MPI_DOUBLE, dest_W, W, lbm_comm);
	}
	if(src_E >= 0)
	{
		MPI_Recv(data, ny, MPI_DOUBLE, src_E, W, lbm_comm, MPI_STATUS_IGNORE);
		for(y=0; y<ny; y++)
			m[y][nx-1] = data[y];
	}
	
	free(data);
}

//shift the matrix m one step to the south
//send the south border to the neighbor
//receive the north border from the neighbor
void shift_S(double **m, int nx, int ny, MPI_Comm lbm_comm)
{
	int y;
	int src_N, dest_S;
	double *data;
	
	MPI_Cart_shift(lbm_comm, 0, 1, &src_N, &dest_S);
	
	data  = malloc(sizeof(double) * (nx));
	
	memcpy(data, &m[ny-1][0], sizeof(double) * (nx));
	
	for(y=(ny-2); y>=0; y--)
	{
		memcpy(&m[y+1][0], &m[y][0], sizeof(double) * (nx));
	}
	
	
	if(dest_S >= 0)
	{
		MPI_Send(data, nx, MPI_DOUBLE, dest_S, S, lbm_comm);
	}
	if(src_N >= 0)
	{
		MPI_Recv(data, nx, MPI_DOUBLE, src_N, S, lbm_comm, MPI_STATUS_IGNORE);
		memcpy(&m[0][0], data, sizeof(double) * (nx));
	}
	
	free(data);
}

//shift the matrix m one step to the north
//send the north border to the neighbor
//receive the south border from the neighbor
void shift_N(double **m, int nx, int ny, MPI_Comm lbm_comm)
{
	int y;
	int src_S, dest_N;
	double *data;
	
	MPI_Cart_shift(lbm_comm, 0, -1, &src_S, &dest_N);
	
	data  = malloc(sizeof(double) * (nx));
	
	memcpy(data, &m[0][0], sizeof(double) * (nx));
	
	for(y=0; y<(ny-1); y++)
	{
		memcpy(&m[y][0], &m[y+1][0], sizeof(double) * (nx));
	}
	
	
	if(dest_N >= 0)
	{
		MPI_Send(data, nx, MPI_DOUBLE, dest_N, N, lbm_comm);
	}
	if(src_S >= 0)
	{
		MPI_Recv(data, nx, MPI_DOUBLE, src_S, N, lbm_comm, MPI_STATUS_IGNORE);
		memcpy(&m[ny-1][0], data, sizeof(double) * (nx));
	}
	
	free(data);
}

void propagation_step(double ***f, int nx, int ny,  MPI_Comm lbm_comm)
{
	//East
	MPI_Barrier(lbm_comm);
	shift_E(f[E], nx, ny, lbm_comm);
	
	//South
	MPI_Barrier(lbm_comm);
	shift_S(f[S], nx, ny, lbm_comm);
	
	//West
	MPI_Barrier(lbm_comm);
	shift_W(f[W], nx, ny, lbm_comm);
	
	//North
	MPI_Barrier(lbm_comm);
	shift_N(f[N], nx, ny, lbm_comm);
	
	//North-East
	MPI_Barrier(lbm_comm);
	shift_NE(f[NE], nx, ny, lbm_comm);
	
	//South-East
	MPI_Barrier(lbm_comm);
	shift_SE(f[SE], nx, ny, lbm_comm);
	
	//South-West
	MPI_Barrier(lbm_comm);
	shift_SW(f[SW], nx, ny, lbm_comm);

	//North-West
	MPI_Barrier(lbm_comm);
	shift_NW(f[NW], nx, ny, lbm_comm);

}
