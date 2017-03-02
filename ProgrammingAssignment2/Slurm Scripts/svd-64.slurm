#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH --nodes=1
# --account=owner-guest
# --partition=kingspeak-guest

#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH --mail-type=END,FAIL 

#SBATCH -o slurm64.out
#SBATCH -e slurm64.err

#module load python/2.7.3
module load matlab

cd $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/

[ ! -d 64-files  ] && mkdir 64-files

cd 64-files

#cp ../matlab-validation.m .
cp ../random-matrix.py .

python random-matrix.py 64 64

srun g++ -o svd.out ../svd.cpp -fopenmp
srun g++ -o svd-seq.out ../svd.cpp

#srun ./svd.out 64 64 -t -p -dm
srun ./svd.out 64 64 -t -n -dm
#matlab -nojvm < matlab-validation.m

#srun ./svd-seq.out 64 64 -t -p -dm
srun ./svd-seq.out 64 64 -t -n -dm
#matlab -nojvm < matlab-validation.m