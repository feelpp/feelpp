% -*- octave -*- 


% Script to generate the Latin Hypercube Sampling for the five inputs
% Kic, D, Q, r and ea

clear all;
clc;

fid = fopen('plhs.dat','w');

n=100;
p=5;

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

D=[Kicmin Dmin Qmin Rmin Eamin; Kicmax Dmax Qmax Rmax Eamax];

AR=lhsdesign(n,p);

for j=1:p
    for i=1:n
        A(i,j) = D(1,j)+AR(i,j)*(D(2,j)-D(1,j));
    end
end

dlmwrite('plhs.dat',A, ' ');

fclose(fid);
