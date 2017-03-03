matrix-size = 128

clean-mat:
	rm *.mat -f

clean: clean-mat
	rm *.out -f
	rm matrix -f

matrix:
	python random-matrix.py $(matrix-size) $(matrix-size)

build: matrix
	g++ -o svd.out svd.cpp -fopenmp
	g++ -o svd-seq.out svd.cpp

run: build
	./svd.out $(matrix-size) $(matrix-size) -t -n -do

runs: build
	./svd-seq.out $(matrix-size) $(matrix-size) -t -n -do 

run-v: build
	./svd.out $(matrix-size) $(matrix-size) -t -p -do

runs-v: build
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