#!/bin/bash

#SBATCH --time=8:00:00
#SBATCH --nodes=1
# --account=owner-guest
# --partition=tangent-guest

# --account=CS4230
# --partition=tangent-freecycle

#SBATCH --account=soc-kp
#SBATCH --partition=soc-kp

#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=kyle.heaton@utah.edu

#SBATCH -o slurm8192-v1.out
#SBATCH -e slurm8192-v1.err

SIZE="8192"
VERSION="-v1"

echo "-------------------------------------------------------------"

cd $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/

[ ! -d $SIZE-files  ] && mkdir $SIZE-files$VERSION
cd $SIZE-files$VERSION

cp ../random-matrix.py .

python random-matrix.py $SIZE $SIZE

srun g++ -o svd$VERSION.out ../svd$VERSION.cpp -fopenmp
srun g++ -o svd-seq.out ../svd-seq.cpp

srun ./svd$VERSION.out $SIZE $SIZE -t -n -dm

srun ./svd-seq.out $SIZE $SIZE -t -n -dm

echo "-------------------------------------------------------------"