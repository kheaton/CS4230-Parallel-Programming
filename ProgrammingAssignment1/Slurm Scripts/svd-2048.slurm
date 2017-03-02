#!/bin/bash

#SBATCH --time=8:00:00
#SBATCH --nodes=1
# --account=owner-guest
# --partition=kingspeak-guest

#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH --mail-type=END,FAIL 

#SBATCH -o slurm2048.out
#SBATCH -e slurm2048.err

#module load python/2.7.3
module load matlab

cd $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/

[ ! -d 2048-files  ] && mkdir 2048-files

cd 2048-files

#cp ../matlab-validation.m .
cp ../random-matrix.py .

python random-matrix.py 2048 2048

srun g++ -o svd.out ../svd.cpp -fopenmp
srun g++ -o svd-seq.out ../svd.cpp

#srun ./svd.out 2048 2048 -t -p -dm
srun ./svd.out 2048 2048 -t -n -dm
#matlab -nojvm < matlab-validation.m

#srun ./svd-seq.out 2048 2048 -t -p -dm
srun ./svd-seq.out 2048 2048 -t -n -dm
#matlab -nojvm < matlab-validation.m