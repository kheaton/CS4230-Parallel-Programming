matrix-size = 512

clean:
	rm *.out -f
	rm matrix -f
	rm *.mat -f

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

run-all: build
	./svd.out 32 32 -t -n -do
	octave -q octave-validation.m
	./svd-seq.out 32 32 -t -n -do
	octave -q octave-validation.m
	./svd.out 64 64 -t -n -do
	octave -q octave-validation.m
	./svd-seq.out 64 64 -t -n -do 
	octave -q octave-validation.m
	./svd.out 128 128 -t -n -do
	octave -q octave-validation.m
	./svd-seq.out 128 128 -t -n -do 
	octave -q octave-validation.m
	./svd.out 256 256 -t -n -do
	octave -q octave-validation.m
	./svd-seq.out 256 256 -t -n -do 
	octave -q octave-validation.m
	./svd.out 512 512 -t -n -do
	octave -q octave-validation.m
	./svd-seq.out 512 512 -t -n -do 
	octave -q octave-validation.m
	./svd.out 1024 1024 -t -n -do
	octave -q octave-validation.m
	./svd-seq.out 1024 1024 -t -n -do 
	octave -q octave-validation.m