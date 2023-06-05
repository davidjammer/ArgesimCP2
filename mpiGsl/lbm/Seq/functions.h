#ifndef __FUNCTIONS
	#define __FUNCTIONS
	
	//Index for 2d matrix
	struct index2d
	{
		int x, y;
	}; typedef struct index2d index2d;
	
	void save(char *file, double **data, int nx, int ny);
	double **create_matrix2d(int x, int y);
	double ***create_matrix3d(int x, int y, int z);
	void free_matrix2d(double **m);
	void free_matrix3d(double ***m, int x, int y, int z);
	void get_driving_ids(double **geometry, int nx, int ny, index2d **driving_ids, int *n_driving_ids);
	void set_value_on_matrix(double **m, int nx, int ny, double value);
	
	int is_number(const char *str);
	int IsPowerOfTwo(int x);
	int set_nx_from_argv(int argc, char **argv, int *Nx);
	int set_ny_from_argv(int argc, char **argv, int *Ny);
	int set_rep_from_argv(int argc, char **argv, int *rep);
	void set_filename_from_argv(int argc, char **argv, char *filename);
	int print_help(int argc, char **argv);
#endif
