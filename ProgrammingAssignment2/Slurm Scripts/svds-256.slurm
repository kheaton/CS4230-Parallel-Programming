#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH -o slurms256-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e slurms256-%j.err-%N # name of the stderr, using job and first node values

# load appropriate modules, in this case Intel compilers, MPICH2
module load openmpi
module load python/2.7.3

# run the program
# see below for ways to do this for different MPI distributions
python $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/random-matrix.py 256 256
mpirun -np $SLURM_NTASKS $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/svd-seq.out 256 256 -t -n -do