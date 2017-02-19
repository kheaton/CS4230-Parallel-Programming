#!/bin/bash

#SBATCH --time=0:5:00
#SBATCH --nodes=1
#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH -o slurm128-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e slurm128-%j.err-%N # name of the stderr, using job and first node values

# load appropriate modules, in this case Intel compilers, MPICH2
module load openmpi
module load python/2.7.3

# run the program
# see below for ways to do this for different MPI distributions
python $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/random-matrix.py 128 128
mpirun -np $SLURM_NTASKS $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/svd.out 128 128 -t -n -do