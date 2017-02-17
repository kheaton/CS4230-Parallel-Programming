#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH -o slurm128-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e slurm128-%j.err-%N # name of the stderr, using job and first node values

module load [gcc, octave]

python random-matrix.py 128 128

./svd.out 128 128 -t -n -do
octave -q octave-validation.m