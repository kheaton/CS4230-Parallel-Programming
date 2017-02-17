#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH -o slurm-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e slurm-%j.err-%N # name of the stderr, using job and first node values

module load [gcc, octave]

python random-matrix.py 512 512

./svd.out 512 512 -t -n -do
octave -q octave-validation.m

./svd-seq.out 512 512 -t -n -do 
octave -q octave-validation.m

python random-matrix.py 4096 4096

./svd.out 4096 4096 -t -n -do
octave -q octave-validation.m

./svd-seq.out 4096 4096 -t -n -do 
octave -q octave-validation.m