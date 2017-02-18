#!/bin/bash

#SBATCH --time=0:5:00
#SBATCH --nodes=1
#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH -o slurm32-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e slurm32-%j.err-%N # name of the stderr, using job and first node values

#module load gcc octave

#python random-matrix.py 32 32

#./svd.out 32 32 -t -n -do
#octave -q octave-validation.m

# set data and working directories
setenv WORKDIR $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2

# load appropriate modules, in this case Intel compilers, MPICH2
module load [gcc,octave,python] openmpi

# run the program
# see below for ways to do this for different MPI distributions
oythong $WORKDIR/random-matrix.py
mpirun -np $SLURM_NTASKS $WORKDIR/svd.out 32 32 -t -n -do
octave -q $WORKDIR/octave-validation.m