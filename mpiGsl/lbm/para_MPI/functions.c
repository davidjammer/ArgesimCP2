#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <mpi.h>

#include <functions.h>
#include <lbm_parameter.h>

int max(int a, int b)
{
	if(a > b)
		return a;
	else
		return b;
}

//save the 2d matrix data with dimesion nx x ny in a file
void save(char *file, double **data, int nx, int ny)
{
	FILE *fd;
	int x, y;
	
	if(access(file, F_OK) == 0)
		remove(file);
	
	fd=fopen(file,"w");
	
	for(y=0; y<ny; y++)
	{
		for(x=0; x<nx; x++)
		{
			fprintf(fd,"%E ",data[y][x]);
		}
		fprintf(fd,"\n");
	}
	fclose(fd);

}

//create a 2d double matrix with dimension x(columns) and y(rows)
double **create_matrix2d(int x, int y)
{
	int i;
	double **m;
	
	m = malloc(sizeof(double*) * y);
	m[0] = calloc(sizeof(double), x * y);
	for(i=0; i<y; i++)
		m[i] = &m[0][i * x];

	return m;
}

void free_matrix2d(double **m)
{
	free(m[0]);	
	free(m);
}

//create a 3d double matrix with dimension x(rolumns), y(rows) and z
double ***create_matrix3d(int x, int y, int z)
{
	int i, j;
	double ***m;
	double *data;
	
	data = calloc(sizeof(double), x * y * z);
	
	m = malloc(sizeof(double**) * z);
	for(j=0; j<z; j++)
	{
		m[j] = malloc(sizeof(double*) * y);
		for(i=0; i<y; i++)
			m[j][i] = &data[i * x + j * x * y];
	}
	return m;
}

void free_matrix3d(double ***m, int x, int y, int z)
{
	int i;

	free(m[0][0]);

	for(i=0; i<z; i++)
		free(m[i]);

	free(m);
}

/*search in the matrix geometry the coordinates of the driving cells
* and put it in to the matrix driving_ids
* Input: geometry -> geometry matrix
*        nx -> size in x direction
*        ny -> size in y direction
* Output: driving_ids -> list of coordinates
*         n_driving_ids -> number of coordinates
*/
void get_driving_ids(double **geometry, int nx, int ny, index2d **driving_ids, int *n_driving_ids)
{
	int x, y;
	
	*n_driving_ids = 0;

	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			if(geometry[y][x] == DRIVING_CELL)
				(*n_driving_ids)++;
	
	*driving_ids = malloc(sizeof(index2d) * (*n_driving_ids));
	
	*n_driving_ids = 0;
	for(y=0; y<ny; y++)
		for(x=0; x<nx; x++)
			if(geometry[y][x] == DRIVING_CELL)
			{
				(*driving_ids)[(*n_driving_ids)].y = y;
				(*driving_ids)[(*n_driving_ids)].x = x;
				(*n_driving_ids)++;
			}

}

//set the value on all elements of matrix m with dimesion nx x ny
void set_value_on_matrix(double **m, int nx, int ny, double value)
{
	int x, y;
	
	for(x=0; x<nx; x++)
		for(y=0; y<ny; y++)
			m[y][x] = value;

}

int is_number(const char *str)
{
	int result = 0;
		
	while(*str !='\0')
	{
		if(*str >= '0' && *str <= '9')
		{
			str++;
			result = 1;
		}
		else
			return 0;
	}
	return result;
}

int IsPowerOfTwo(int x)
{
    return (!(x & (x - 1)) && x);
}

//set parameter Nx from argv
int set_nx_from_argv(int argc, char **argv, int *Nx)
{
	int number;
	int i = 0;
	
	do
	{
		if(strcmp(argv[i], "-Nx") == 0)
		{
			if(is_number(argv[i+1]))
			{
				number = atoi(argv[i+1]);
				if(IsPowerOfTwo(number))
				{
					*Nx = number;
					return 0;
				}
				else
					return -1;
			}
			else
				return -1;
		}
	}while(argv[++i] != NULL);

	return 0;
}

//set parameter Ny from argv
int set_ny_from_argv(int argc, char **argv, int *Ny)
{
	int number;
	int i = 0;
	
	do
	{
		if(strcmp(argv[i], "-Ny") == 0)
		{
			if(is_number(argv[i+1]))
			{
				number = atoi(argv[i+1]);
				if(IsPowerOfTwo(number))
				{
					*Ny = number;
					return 0;
				}
				else
					return -1;
			}
			else
				return -1;
		}
	}while(argv[++i] != NULL);
	
	return 0;
}

//set parameter Rep from argv
int set_rep_from_argv(int argc, char **argv, int *rep)
{
	int number;
	int i = 0;
	
	do
	{
		if(strcmp(argv[i], "-Rep") == 0)
		{
			if(is_number(argv[i+1]))
			{
				number = atoi(argv[i+1]);
				
				*rep = number;
				return 0;
			}
			else
				return -1;
		}
	}while(argv[++i] != NULL);

	return 0;
}

void set_filename_from_argv(int argc, char **argv, char *filename)
{
	int i = 0;
	
	do
	{
		if(strcmp(argv[i], "-o") == 0)
		{
			if(argv[i+1] != NULL)
				sprintf(filename, "%s", argv[i+1]);
				return;
		}
	}while(argv[++i] != NULL);

}

int print_help(int argc, char **argv)
{
	int world_rank;
	
	int i = 0;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	do
	{
		if(strcmp(argv[i], "--help") == 0)
		{
			if(world_rank == 0)
			{
				printf("usage: mpirun -np <Number of processes> ./lbm [OPTIONS]\n");
				printf("Options:\n");
				printf("-Nx\t\tNumber of elelments in x direction (default: 256)\n");
				printf("-Ny\t\tNumber of elelments in y direction (default: 256)\n");
				printf("-Rep\t\tNumber of repetitions (default: 350000)\n"); 
				printf("-o\t\tOutput file name (default: u.dat))\n");
			}
			return 1;
		}
	}while(argv[++i] != NULL);
	return 0;
}
