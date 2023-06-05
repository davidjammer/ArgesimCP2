# ARGESIM Benchmark on Parallel and Distributed Simulation / Case Study 3 – Solution of a Partial Differential Equation
Solution of task C3 of the SNE benchmark CP2. Definition: https://www.sne-journal.org/fileadmin/user_upload_sne/benchmarks/CP/CP2-definition.pdf

* Dependencies:
  * libgsl-dev
  * openmpi

* Build: make all

* parallel implementation:
  * usage: mpirun -np (number of processes) ./para [N] [dt] [STEPS]
  * N = number of intervals
  * dt = step size
  * STEPS = number of steps
  * A call without passing parameters starts a simulation with the default parameters <br>
  (N = 500, dt = 0.001, STEPS = 10000)
* sequential implemetation:
  * usage: ./seq [N] [dt] [STEPS]
  * N = number of intervals
  * dt = step size
  * STEPS = number of steps
  * A call without passing parameters starts a simulation with the default parameters <br>
  (N = 500, dt = 0.001, STEPS = 10000)
* output: The sequential and parallel versions create a file named "u.dat". <br>
The first row contains the location -> u.dat(1, 2:end) <br>
The first column contains the time -> u.dat(:, 2:end) <br>
The excitation over the space is one row per time and results from the rest of the file -> u.dat(2:end, 2:end)
* plot:
  * plotResults("u.dat")
