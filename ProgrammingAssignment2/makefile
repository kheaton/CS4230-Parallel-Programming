matrix-size = 128
matrix-size2 = 512

clean:
	rm *.out -f
	rm matrix -f
	rm *.mat -f

matrix:
	python random-matrix.py $(matrix-size) $(matrix-size)

build: matrix
	g++ -o svd.out svd.cpp -fopenmp

run: build
	./svd.out $(matrix-size) $(matrix-size) -t -n -do

run-v: build
	./svd.out $(matrix-size) $(matrix-size) -t -p -do

rebuild: clean build

valid: matrixA.mat matrixScpu.mat matrixUcpu.mat matrixVcpu.mat
	octave -q octave-validation.m