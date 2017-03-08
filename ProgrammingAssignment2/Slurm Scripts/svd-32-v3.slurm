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

#SBATCH -o slurm32-v2.out
#SBATCH -e slurm32-v2.err

SIZE="32"
VERSION="-v3"

echo "-------------------------------------------------------------"

cd $HOME/CS4230-Parallel-Programming/ProgrammingAssignment2/

[ ! -d ./Files/$SIZE-files$VERSION ] && mkdir ./Files/$SIZE-files$VERSION
cd ./Files/$SIZE-files$VERSION

cp ../../random-matrix.py .

python random-matrix.py $SIZE $SIZE

srun g++ -o svd$VERSION.out ../../svd$VERSION.cpp -fopenmp
srun g++ -o svd-seq.out ../../svd-seq.cpp

srun ./svd$VERSION.out $SIZE $SIZE -t -n -dm -1
srun ./svd$VERSION.out $SIZE $SIZE -t -n -dm 2
srun ./svd$VERSION.out $SIZE $SIZE -t -n -dm 4
srun ./svd$VERSION.out $SIZE $SIZE -t -n -dm 8
srun ./svd$VERSION.out $SIZE $SIZE -t -n -dm 16
srun ./svd$VERSION.out $SIZE $SIZE -t -n -dm 32
srun ./svd$VERSION.out $SIZE $SIZE -t -n -dm 64
srun ./svd$VERSION.out $SIZE $SIZE -t -n -dm 128

srun ./svd-seq.out $SIZE $SIZE -t -n -dm

echo "-------------------------------------------------------------"