#ifndef __LBM_PROPAGATION_STEP
	#define __LBM_PROPAGATION_STEP
	
	void propagation_step(double ***f, int nx, int ny);
	void shift_N(double **m, int nx, int ny);
	void shift_S(double **m, int nx, int ny);
	void shift_W(double **m, int nx, int ny);
	void shift_E(double **m, int nx, int ny);
	void shift_NE(double **m, int nx, int ny);
	void shift_SE(double **m, int nx, int ny);
	void shift_SW(double **m, int nx, int ny);
	void shift_NW(double **m, int nx, int ny);
#endif
