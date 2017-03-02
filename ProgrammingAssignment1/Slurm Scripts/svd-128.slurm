#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH --nodes=1
# --account=owner-guest
# --partition=kingspeak-guest

#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH --mail-type=END,FAIL

#SBATCH -o slurm128.out
#SBATCH -e slurm128.err

#module load python/2.7.3
module load matlab

cd $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/

[ ! -d 128-files  ] && mkdir 128-files

cd 128-files

#cp ../matlab-validation.m .
cp ../random-matrix.py .

python random-matrix.py 128 128

srun g++ -o svd.out ../svd.cpp -fopenmp
srun g++ -o svd-seq.out ../svd.cpp

#srun ./svd.out 128 128 -t -p -dm
srun ./svd.out 128 128 -t -n -dm
#matlab -nojvm < matlab-validation.m

#srun ./svd-seq.out 128 128 -t -p -dm
srun ./svd-seq.out 128 128 -t -n -dm
#matlab -nojvm < matlab-validation.m