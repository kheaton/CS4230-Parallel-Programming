#!/bin/bash

#SBATCH --time=8:00:00
#SBATCH --nodes=1
# --account=owner-guest
# --partition=kingspeak-guest

#SBATCH --account=CS4230
#SBATCH --partition=kingspeak-freecycle

#SBATCH --mail-type=END,FAIL 

#SBATCH -o slurm4096.out
#SBATCH -e slurm4096.err

#module load python/2.7.3
module load matlab

cd $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/

[ ! -d 4096-files  ] && mkdir 4096-files

cd 4096-files

#cp ../matlab-validation.m .
cp ../random-matrix.py .

python random-matrix.py 4096 4096

srun g++ -o svd.out ../svd.cpp -fopenmp
srun g++ -o svd-seq.out ../svd.cpp

#srun ./svd.out 4096 4096 -t -p -dm
srun ./svd.out 4096 4096 -t -n -dm
#matlab -nojvm < matlab-validation.m

#srun ./svd-seq.out 4096 4096 -t -p -dm
srun ./svd-seq.out 4096 4096 -t -n -dm
#matlab -nojvm < matlab-validation.m