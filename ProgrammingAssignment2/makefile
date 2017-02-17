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

buildi: matrix
	icpc -qopenmp -o svd-i.out svd.cpp
	icpc -o svd-seq-i.out svd.cpp

runi: buildi
	./svdi.out $(matrix-size) $(matrix-size) -t -n -do

runsi: buildi
	./svd-seqi.out $(matrix-size) $(matrix-size) -t -n -do 

runi-v: buildi
	./svdi.out $(matrix-size) $(matrix-size) -t -p -do

runsi-v: buildi
	./svd-seqi.out $(matrix-size) $(matrix-size) -t -p -do

run-bothi: runi runsi

run-bothi-v: runi-v runsi-v

rebuildi: clean buildi
