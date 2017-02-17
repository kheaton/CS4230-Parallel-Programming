#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH -o slurm1024-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e slurm1024-%j.err-%N # name of the stderr, using job and first node values

module load [gcc, octave]

python random-matrix.py 1024 1024

./svd.out 1024 1024 -t -n -do
octave -q octave-validation.m