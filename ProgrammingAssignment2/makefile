matrix-size = 128
matrix-size2 = 512

clean:
	rm *.out -f
	rm matrix -f
	rm *.dat -f

matrix:
	python random-matrix.py $(matrix-size) $(matrix-size)

matrix2:
	pythong random-matrix.py $(matrix-size) $(matrix-size2)

build: matrix
	g++ -o svd.out svd.cpp -fopenmp

build2: matrix2
	g++ -o svd.out svd.cpp -fopenmp

run: build
	./svd.out $(matrix-size) $(matrix-size) -t -do

run2: build2
	./svd.out $(matrix-size2) $(matrix-size2) -t -do

run-v: build
	./svd.out $(matrix-size) $(matrix-size) -t -p -do

run-v2: build2
	./svd.out $(matrix-size2) $(matrix-size2) -t -p -do

rebuild: clean build

rebuild2: clean build2

svd-omp: OMP-SVD.cpp
	g++ -o svd-omp.out OMP-SVD.cpp -fopenmp

do-omp: svd-omp
	./svd-omp.out $(matrix-size) $(matrix-size) -t -p -do