#!/bin/bash

#SBATCH --time=8:00:00
#SBATCH --nodes=1
# --account=owner-guest
# --partition=kingspeak-guest

#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH --mail-type=END,FAIL 

#SBATCH -o slurm1024.out
#SBATCH -e slurm1024.err

#module load python/2.7.3
module load matlab

cd $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/

[ ! -d 1024-files  ] && mkdir 1024-files

cd 1024-files

#cp ../matlab-validation.m .
cp ../random-matrix.py .

python random-matrix.py 1024 1024

srun g++ -o svd.out ../svd.cpp -fopenmp
srun g++ -o svd-seq.out ../svd.cpp

#srun ./svd.out 1024 1024 -t -p -dm
srun ./svd.out 1024 1024 -t -n -dm
#matlab -nojvm < matlab-validation.m

#srun ./svd-seq.out 1024 1024 -t -p -dm
srun ./svd-seq.out 1024 1024 -t -n -dm
#matlab -nojvm < matlab-validation.m