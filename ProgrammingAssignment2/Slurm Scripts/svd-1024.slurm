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

#SBATCH -o slurm1024.out
#SBATCH -e slurm1024.err

SIZE="1024"
VERSION=""

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