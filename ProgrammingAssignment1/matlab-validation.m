Ucpu = dlmread('matrixUcpu.dat'); 
Vcpu = dlmread('matrixVcpu.dat'); 
Scpu = dlmread('matrixScpu.dat'); 
A = dlmread('matrixA.dat');

[U,S,V] = svd(A);
ResU = abs(mean(mean(U))); 
ResS = abs(mean(mean(S))); 
ResV = abs(mean(mean(V)));
ResU1 = abs(mean(mean(Ucpu))); 
ResS1 = abs(mean(mean(Scpu))); 
ResV1 = abs(mean(mean(Vcpu)));
ru = abs(ResU - ResU1) rv = abs(ResV - ResV1) rs = abs(ResS - ResS1)