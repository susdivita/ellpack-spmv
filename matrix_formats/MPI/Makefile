MPI_CC=mpic++
CC=g++
CFLAGS=-lstdc++ -ldl -lm -llsb
CFLAGS_EULER=-lstdc++ -ldl -lm -fopenmp
EXE=./smvm
SRC=main.cpp DDDNaiveDistributed.cpp

build:
	$(MPI_CC) $(SRC) $(CFLAGS) -O3 -o $(EXE)

debug-build:
	$(MPI_CC) $(SRC) $(CFLAGS) -O0 -g -o $(EXE)_g

build_euler:
	$(MPI_CC) -std=c++11 $(SRC) $(CFLAGS_EULER) -O3 -I ~/liblsb-0.2.2/build/include -L ~/liblsb-0.2.2/build/lib -llsb -o $(EXE)

# Usage: make run PROCS=2 FORMAT=DDDNaive \
# MTX_IN="../../inputs/Materials\ Problem/arc130.mtx" \
# VEC_IN="../../inputs/Materials\ Problem/test_vec_130_ones_in.txt"
run: $(EXE)
	mpirun --np $(PROCS) $(EXE) $(FORMAT) $(MTX_IN) $(VEC_IN)

# Usage: srun --cpus-per-task=N make run_euler FORMAT=. MTX_IN=.. VEC_IN=..
run_euler:
	$(EXE) $(FORMAT) $(MTX_IN) $(VEC_IN)

debug-valgrind:
	valgrind --leak-check=yes mpirun --np 1 $(EXE) $(FORMAT) $(MTX_IN) $(VEC_IN)

debug-gdb:
	mpirun --np 1 gdb $(EXE) $(FORMAT) $(MTX_IN) $(VEC_IN)
