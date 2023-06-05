para_MPI: parallele Lösung mit MPI
mpicc -O3 main.c -lm
mpirun -np 4 lbm -Nx 256 -Ny 256 -Rep 10000 -o u.dat

Seq: seq. Lösung
gcc -O3 main.c -lm


