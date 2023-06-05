INC = .
OBJ = functions.o lbm_init.o lbm_collision_step.o lbm_propagation_step.o lbm_magnitude.o 
MAINOBJ = main.o test_shift.o
CC = mpicc
CFLAGS = -O3
LDFLAGS = -lm

all: lbm 

lbm: $(OBJ) $(MAINOBJ) 
	@echo [CC] lbm
	@$(CC) $(CFLAGS) -I $(INC) $(OBJ) main.o $(LDFLAGS) -o lbm
	@echo [CC] test_shit
	@$(CC) $(CFLAGS) -I $(INC) $(OBJ) test_shift.o $(LDFLAGS) -o test_shift

%.o: %.c
	@echo [CC] $<
	@$(CC) $(CFLAGS) -I $(INC) $(LDFLAGS) -c $< -o $@

clean:
	rm -rf lbm test_shift *.dat *.o *~
