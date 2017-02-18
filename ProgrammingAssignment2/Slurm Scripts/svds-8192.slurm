#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH -o slurm8192-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e slurm8192-%j.err-%N # name of the stderr, using job and first node values

module load [gcc, octave]

python random-matrix.py 8192 8192

./svd-seq.out 8192 8192 -t -n -do
octave -q octave-validation.m