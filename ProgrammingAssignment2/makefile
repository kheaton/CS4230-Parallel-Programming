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