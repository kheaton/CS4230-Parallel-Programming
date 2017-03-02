#!/bin/bash

#SBATCH --time=8:00:00
#SBATCH --nodes=1
# --account=owner-guest
# --partition=kingspeak-guest

#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH --mail-type=END,FAIL 

#SBATCH -o slurm8192.out
#SBATCH -e slurm8192.err

#module load python/2.7.3
module load matlab

cd $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/

[ ! -d 8192-files  ] && mkdir 8192-files

cd 8192-files

#cp ../matlab-validation.m .
cp ../random-matrix.py .

python random-matrix.py 8192 8192

srun g++ -o svd.out ../svd.cpp -fopenmp
srun g++ -o svd-seq.out ../svd.cpp

#srun ./svd.out 8192 8192 -t -p -dm
srun ./svd.out 8192 8192 -t -n -dm
#matlab -nojvm < matlab-validation.m

#srun ./svd-seq.out 8192 8192 -t -p -dm
srun ./svd-seq.out 8192 8192 -t -n -dm
#matlab -nojvm < matlab-validation.m