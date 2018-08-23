%  -*-matlab-*-

% Script for the SA using Latin Hypercube Sampling ...

clear all;
clc;

fid = fopen('mplhs.dat','w');

n=300;
p=3;
% Input variation domains---
Kicmin = 2.0e-1;
Kicmax = 1.5e+2;

Dmin   = 5.0e-4;
Dmax   = 1.0e-2;

Qmin   = 0;
Qmax   = 1.0e+6;

Rmin   = 1.0e-1;
Rmax   = 1.0e+2;

Eamin  = 2.5e-3;
Eamax  = 5.0e-2;
%---------------------------

D=[Kicmin Dmin Rmin; Kicmax Dmax Rmax];

AR=lhsdesign(n,p);

for j=1:p
    for i=1:n
        P(i,j) = D(1,j)+AR(i,j)*(D(2,j)-D(1,j));
    end
end

dlmwrite('mplhs.dat', P, ' ');


fclose(fid);
