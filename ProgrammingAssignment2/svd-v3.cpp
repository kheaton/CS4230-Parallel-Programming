// svd.cpp with SIMD optimizations

#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sys/time.h>
#include "omp.h"
#include<cstdlib>

#define epsilon 1.e-8

using namespace std;

template <typename T> double sgn(T val)
{
	return (val > T(0)) - (val < T(0));
}

int main (int argc, char* argv[]) {

	int M,N;
	bool octave; 
	string T,P,Db, threads;
	M = atoi(argv[1]);
	N = atoi(argv[2]);

	double elapsedTime,elapsedTime2;
	timeval start,end,end2;

	if(argc < 4) {
		cout<<"Please input the size of Matrix and at least one of the options: -t -p -do/-dm";
		return 0;
	}

	if(M != N) {
		cout<<"Error: Matrix must be square";
		return 0;
	}
  
	if(argc > 3) {
    	T = argv[3];
    	if(argc > 4) {
      		P = argv[4];
      		if(argc > 5) {
        		Db = argv[5];
				if(Db=="-do")
				{
					cout << "Using Octave" << endl;
					octave = true;
				}
				else if(Db=="-dm")
				{
					cout << "Using Matlab" << endl;
					octave = false;
				}
				else{
					cout << "Invalid debug option %s" << Db <<  endl;
					return 1;
				}

				threads = argv[6];
				cout << "Threads :" << threads << endl;
			}
    	}
  	}
  
	int numberOfThreads = std::atoi (threads.c_str());

	if(numberOfThreads > 0) {
		omp_set_dynamic(0);
		omp_set_num_threads(numberOfThreads);
	}
	
	double **U,**V, *S,**U_t, **V_t, **A;
  	double alpha, beta, gamma, c, zeta, t,s,sub_zeta, converge;
	double t1, t2, t3, t4;

  	int acum = 0;
  	int temp1, temp2;
  	converge = 1.0;

  	U = new double*[N];
  	V = new double*[N];
  	U_t = new double*[N];
  	V_t = new double*[N];
  	A = new double*[N];
  	S = new double[N];

  	for(int i =0; i<N; i++) {
		U[i] = new double[N];
 		V[i] = new double[N];
		U_t[i] = new double[N];
		V_t[i] = new double[N];
		A[i] = new double[N];
  	}

  	//Read from file matrix, if not available, app quit
  	//Already transposed

  	ifstream matrixfile("matrix");
  	if(!(matrixfile.is_open())) {
    	cout<<"Error: file not found"<<endl;
    	return 0;
  	}

  	for(int i = 0; i < M; i++) {
    	for(int j =0; j < N; j++) {
      		matrixfile >> U_t[i][j];
    	}
  	}

  	matrixfile.close();
 
  	for(int i=0; i<M;i++) {
    	for(int j=0; j<N;j++) {
			if(i==j) {
				V_t[i][j] = 1.0;
			}
			else {
				V_t[i][j] = 0.0;
			}
		}
  	}

  	//Store A for debug purpouse
   	for(int i=0; i<M;i++) {
      	for(int j=0; j<N;j++) {
       		A[i][j] = U_t[j][i];
      	}
	}

  	gettimeofday(&start, NULL);

   	double conv;
   	while(converge > epsilon) { 		//convergence
    	converge = 0.0;	
    	acum++;							//counter of loops

    	for(int i = 1; i < M; i++) {
      		for(int j = 0; j < i; j++) {

          		alpha = 0.0;
          		beta = 0.0;
          		gamma = 0.0;

          		for(int k = 0; k < N; k++) {
            		alpha = alpha + (U_t[i][k] * U_t[i][k]);
            		beta = beta + (U_t[j][k] * U_t[j][k]);
            		gamma = gamma + (U_t[i][k] * U_t[j][k]);
          		}

          		converge = max(converge, abs(gamma) / sqrt(alpha * beta));	//compute convergence
	  																		//basicaly is the angle
																			//between column i and j


				zeta = (beta - alpha) / (2.0 * gamma);
				t = sgn(zeta) / (abs(zeta) + sqrt(1.0 + (zeta*zeta)));	//compute tan of angle
				c = 1.0 / (sqrt (1.0 + (t * t)));							//extract cos
				s = c * t;												//extrac sin
 

	  			//Apply rotations on U and V
				for(int k=0; k<N; k += 4) {
					
					t = U_t[i][k];
					U_t[i][k] = c * t - s * U_t[j][k];
					U_t[j][k] = s * t + c * U_t[j][k];

					t = V_t[i][k];
					V_t[i][k] = c * t - s * V_t[j][k];
					V_t[k][k] = s * t + c * V_t[j][k];

					t1 = U_t[i][k+1];
					U_t[i][k+1] = c * t1 - s * U_t[j][k+1];
					U_t[j][k+1] = s * t1 + c * U_t[j][k+1];

					t1 = V_t[i][k+1];
					V_t[i][k+1] = c * t1 - s * V_t[j][k+1];
					V_t[k][k+1] = s * t1 + c * V_t[j][k+1];

					t2 = U_t[i][k+2];
					U_t[i][k+2] = c * t2 - s * U_t[j][k+2];
					U_t[j][k+2] = s * t2 + c * U_t[j][k+2];

					t2 = V_t[i][k+2];
					V_t[i][k+2] = c * t2 - s * V_t[j][k+2];
					V_t[k][k+2] = s * t2 + c * V_t[j][k+2];

					t3 = U_t[i][k+3];
					U_t[i][k+3] = c * t3 - s * U_t[j][k+3];
					U_t[j][k+3] = s * t3 + c * U_t[j][k+3];

					t1 = V_t[i][k+3];
					V_t[i][k+3] = c * t3 - s * V_t[j][k+3];
					V_t[k][k+3] = s * t3 + c * V_t[j][k+3];

					/*
					t = U_t[i][k];
					U_t[i][k] = c * t - s * U_t[j][k];
					U_t[j][k] = s * t + c * U_t[j][k];

					t = V_t[i][k];
					V_t[i][k] = c * t - s * V_t[j][k];
					V_t[j][k] = s * t + c * V_t[j][k];
					*/
				}
			}
		}
	}

  	//Create matrix S
	#pragma openmp parallel for private(t)
	for(int i =0; i<M; i++) {
		t=0;
    	for(int j=0; j<N;j++) {
      		t=t + pow(U_t[i][j],2);
    	}
    	t = sqrt(t);

		#pragma openmp parallel for
		for(int j=0; j<N;j += 4) {
      		U_t[i][j] = U_t[i][j] / t;
      		U_t[i][j+1] = U_t[i][j+1] / t;
      		U_t[i][j+2] = U_t[i][j+2] / t;
      		U_t[i][j+3] = U_t[i][j+3] / t;
			  
      		if(i == j) {
        		S[i] = t;
      		}
    	}
  	}

  	gettimeofday(&end, NULL);

 	/************************************************************/
  	for(int i =0; i<M; i++) {
    	for(int j =0; j<N; j++) {
      		U[i][j] = U_t[j][i];
      		V[i][j] = V_t[j][i];
		}
  	}

	//Output time and iterations
  	if(T=="-t" || P =="-t") {
		cout << "Parallel: " << numberOfThreads << endl;
    	cout << "iterations: " << acum << endl;
    	elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
    	elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
    	cout << "Time: " << elapsedTime << " ms." << endl << endl;
  	}

  	// Output the matrixes for debug
  	if(T== "-p" || P == "-p") {
  		cout << "U" << endl << endl;
  		for(int i =0; i<M; i++) {
    		for(int j =0; j<N; j++) {
      			cout<<U[i][j]<<"  ";
    		}
    		cout<<endl;
  		}

  		cout << endl << "V" << endl << endl;
  		for(int i =0; i<M; i++) {
    		for(int j =0; j<N; j++) {
      			cout << V[i][j] << "  ";
    		}
    		cout << endl;
  		}

  		cout << endl << "S" << endl << endl;
  		for(int i =0; i<M; i++) {
    		for(int j =0; j<N; j++) {
       			if(i==j){ cout << S[i] << "  "; }
       			else { cout << "0.0  "; }
    		}
			cout << endl;
  		}
  	}

	//Generate Octave files for debug purpouse
   	if(Db == "-do" || T == "-do" || P == "-do" || Db == "-dm" || T == "-dm" || P == "-dm") {
    	ofstream Af;
    	//file for Matrix A
		if(octave) {
			Af.open("matrixA.mat"); 
			Af << "# Created from debug\n# name: A\n# type: matrix\n# rows: " << M << "\n# columns: " << N <<"\n";
		}
		else {
			Af.open("matrixA.dat");
		}

    	for(int i = 0; i<M;i++) {
      		for(int j =0; j<N;j++) {
        		Af << " " << A[i][j];
      		}
      		Af<<"\n";
    	}
    
    	Af.close();

    	ofstream Uf;

    	//File for Matrix U
		if(octave) {
			Uf.open("matrixUcpu.mat");
			Uf << "# Created from debug\n# name: Ucpu\n# type: matrix\n# rows: " << M <<"\n# columns: " << N <<"\n";
		}
		else {
			Uf.open("matrixUcpu.dat");
		}
    	for(int i = 0; i<M;i++) {
      		for(int j =0; j<N;j++) {
				Uf << " " << U[i][j];
      		}
      		Uf<<"\n";
    	}
    	Uf.close();

    	ofstream Vf;
    	//File for Matrix V
		if(octave) {
	    	Vf.open("matrixVcpu.mat");
			Vf << "# Created from debug\n# name: Vcpu\n# type: matrix\n# rows: " << M << "\n# columns: " << N <<"\n";
		}
		else {
			Vf.open("matrixVcpu.dat");
		}

    	for(int i = 0; i<M;i++) {
    		for(int j =0; j<N;j++) {
        		Vf << " " << V[i][j];
      		}
			Vf<<"\n";
    	}
    
    	Vf.close();

    	ofstream Sf;
    	//File for Matrix S
    	if(octave) {
			Sf.open("matrixScpu.mat");
			Sf << "# Created from debug\n# name: Scpu\n# type: matrix\n# rows: " << M << "\n# columns: " << N <<"\n";
		}
		else {
			Sf.open("matrixScpu.dat");
		}
    
    	for(int i = 0; i<M;i++) {
      		for(int j =0; j<N;j++) {
        		if(i == j) {
         			Sf << " " << S[i];
        		}
        		else {
          			Sf << " 0.0";
        		}
      		}
      		Sf << "\n";
    	}
    
    	Sf.close();
 	}

   	delete [] S;
   	for(int i = 0; i<N;i++) {
		delete [] A[i];
	   	delete [] U[i];
	   	delete [] V[i];
	   	delete [] U_t[i];
	   	delete [] V_t[i];
	}

	return 0;
}