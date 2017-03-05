matrix-size = 512

clean-err-results:
	rm slurm*.err

clean-out-results:
	rm slurm*.out

clean-results: clean-out-results clean-err-results

clean-mat:
	rm *.mat -f

clean: clean-mat
	rm *^slurm*.out -f
	rm matrix -f

matrix:
	python random-matrix.py $(matrix-size) $(matrix-size)

build: matrix
	g++ -o svd.out svd.cpp -fopenmp

build-seq: matrix
	g++ -o svd-seq.out svd-seq.cpp

build-all: matrix build-seq build
	
run: build
	./svd.out $(matrix-size) $(matrix-size) -t -n -do

runs: build-seq
	./svd-seq.out $(matrix-size) $(matrix-size) -t -n -do 

run-v: build-all
	./svd.out $(matrix-size) $(matrix-size) -t -p -do

runs-v: build-all
	./svd-seq.out $(matrix-size) $(matrix-size) -t -p -do

run-both: run runs

run-both-v: run-v runs-v

rebuild: clean build

valid: matrixA.mat matrixScpu.mat matrixUcpu.mat matrixVcpu.mat
	octave -q octave-validation.m

validm: matrixA.mat matrixScpu.mat matrixUcpu.mat matrixVcpu.mat
	matlab -nojvm < matlab-validation.m

validate: run valid clean-mat

validates: runs valid clean-mat

validate-both: validate validates

build1: matrix
	g++ -o svd-v1.out svd-v1.cpp -fopenmp

run1: build1
	./svd-v1.out $(matrix-size) $(matrix-size) -t -n -do

build2: matrix
	g++ -o svd-v2.out svd-v2.cpp -fopenmp

run2: build2
	./svd-v2.out $(matrix-size) $(matrix-size) -t -n -do

build3: matrix
	g++ -o svd-v3.out svd-v3.cpp -fopenmp

run3: build3
	./svd-v3.out $(matrix-size) $(matrix-size) -t -n -do

build4: matrix
	g++ -o svd-v4.out svd-v4.cpp -fopenmp

run4: build4
	./svd-v4.out $(matrix-size) $(matrix-size) -t -n -do

batch-all:
	sbatch  ./Slurm\ Scripts/svd-32.slurm
	sbatch  ./Slurm\ Scripts/svd-32-v1.slurm
	sbatch  ./Slurm\ Scripts/svd-32-v2.slurm
	sbatch  ./Slurm\ Scripts/svd-64.slurm
	sbatch  ./Slurm\ Scripts/svd-64-v1.slurm
	sbatch  ./Slurm\ Scripts/svd-64-v2.slurm
	sbatch  ./Slurm\ Scripts/svd-128.slurm
	sbatch  ./Slurm\ Scripts/svd-128-v1.slurm
	sbatch  ./Slurm\ Scripts/svd-128-v2.slurm
	sbatch  ./Slurm\ Scripts/svd-256.slurm
	sbatch  ./Slurm\ Scripts/svd-256-v1.slurm
	sbatch  ./Slurm\ Scripts/svd-256-v2.slurm
	sbatch  ./Slurm\ Scripts/svd-512.slurm
	sbatch  ./Slurm\ Scripts/svd-512-v1.slurm
	sbatch  ./Slurm\ Scripts/svd-512-v2.slurm
	sbatch  ./Slurm\ Scripts/svd-1024.slurm
	sbatch  ./Slurm\ Scripts/svd-1024-v1.slurm
	sbatch  ./Slurm\ Scripts/svd-1024-v2.slurm
	sbatch  ./Slurm\ Scripts/svd-2048.slurm
	sbatch  ./Slurm\ Scripts/svd-2048-v1.slurm
	sbatch  ./Slurm\ Scripts/svd-2048-v2.slurm
	sbatch  ./Slurm\ Scripts/svd-4096.slurm
	sbatch  ./Slurm\ Scripts/svd-4096-v1.slurm
	sbatch  ./Slurm\ Scripts/svd-4096-v2.slurm
	sbatch  ./Slurm\ Scripts/svd-8192.slurm
	sbatch  ./Slurm\ Scripts/svd-8192-v1.slurm
	sbatch  ./Slurm\ Scripts/svd-8192-v2.slurm