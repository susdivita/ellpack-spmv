#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=6
#SBATCH --ntasks-per-socket=2
#SBATCH --cores-per-socket=12
#SBATCH --nodes=2
export OMP_NUM_THREADS=6
export OMP_PLACES=cores
srun --cpus-per-task=6 make run_euler FORMAT=DDDNaive MTX_IN="../../inputs/Materials\ Problem/arc130.mtx" VEC_IN="../../inputs/Materials\ Problem/test_vec_130_random_in.txt"

