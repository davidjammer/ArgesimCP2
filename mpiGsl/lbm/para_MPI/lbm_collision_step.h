#ifndef __LBM_COLLISION_STEP
	#define __LBM_COLLISION_STEP
	
	void calc_rho(double **rho, double ***f, int nx, int ny);
	void calc_ux(double **ux, double ***f, double **rho, int nx, int ny);
	void calc_uy(double **uy, double ***f, double **rho, int nx, int ny);
	void set_ux_for_driving_cells(double **ux, index2d *driving_ids, int n_driving_ids);
	void set_uy_for_driving_cells(double **uy, index2d *driving_ids, int n_driving_ids);
	void calc_usqr(double **usqr, double **ux, double **uy, int nx, int ny);
	void calc_feq(double ***feq, double **rho, double **ux, double **uy, double **usqr, int nx, int ny);
	void f_cell_calc(double ***f, double ***feq, double **geometry, double tau, int nx, int ny);
#endif
