load matrixUcpu.mat;
load matrixVcpu.mat;
load matrixScpu.mat;
load matrixA;

[U,S,V] = svd(A);

ResU = abs(mean(mean(U)));
ResS = abs(mean(mean(S)));
ResV = abs(mean(mean(V)));

ResU1 = abs(mean(mean(Ucpu)));
ResS1 = abs(mean(mean(Scpu)));
ResV1 = abs(mean(mean(Vcpu)));

ru = abs(ResU - ResU1)
rv = abs(ResV - ResV1)
rs = abs(ResS - ResS1)

