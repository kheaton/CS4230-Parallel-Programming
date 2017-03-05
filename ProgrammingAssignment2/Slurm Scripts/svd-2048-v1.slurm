#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH --nodes=1
# --account=owner-guest
# --partition=tangent-guest

# --account=CS4230
# --partition=tangent-freecycle

#SBATCH --account=soc-kp
#SBATCH --partition=soc-kp

#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=kyle.heaton@utah.edu

#SBATCH -o slurm2048-v1.out
#SBATCH -e slurm2048-v1.err

SIZE="2048"
VERSION="-v1"

cd $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/

[ ! -d $SIZE-files  ] && mkdir $SIZE-files$VERSION
cd $SIZE-files$VERSION

#cp ../matlab-validation.m .
cp ../random-matrix.py .

python random-matrix.py $SIZE $SIZE

srun g++ -o svd$VERSION.out ../svd$VERSION.cpp -fopenmp
srun g++ -o svd-seq.out ../svd-seq.cpp

srun ./svd$VERSION.out $SIZE $SIZE -t -n -dm
#matlab -nojvm < matlab-validation.m

srun ./svd-seq.out $SIZE $SIZE -t -n -dm
#matlab -nojvm < matlab-validation.m