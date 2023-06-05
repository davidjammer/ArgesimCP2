INC = .
OBJ = functions.o lbm_init.o lbm_collision_step.o lbm_propagation_step.o lbm_magnitude.o main.o
CC = gcc
CFLAGS = -O3
LDFLAGS = -lm

all: lbm

lbm: $(OBJ)
	@echo [CC] lbm
	@$(CC) $(CFLAGS) -I $(INC) $(OBJ) $(LDFLAGS) -o lbm
%.o: %.c
	@echo [CC] $<
	@$(CC) $(CFLAGS) -I $(INC) $(LDFLAGS) -c $< -o $@

clean:
	rm -rf lbm *.dat *.o
