#ifndef __PARAMETER
	#define __PARAMETER
	
	#define RHO_0 1.0
	#define u_0 0.1
	#define RE 1000
	
	
	#define DEFAULT_Nx 256
	#define DEFAULT_Ny 256
	#define DEFAULT_REP 350000
	
	//Constants for direction
	enum direction
	{
		C=0,
		E,
		S,
		W,
		N,
		NE,
		SE,
		SW,
		NW,
	};
	
	#define FLUID_CELL 0.0
	#define WALL_CELL 1.0
	#define DRIVING_CELL 2.0
	
#endif
