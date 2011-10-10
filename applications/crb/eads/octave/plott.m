# -*- octave -*- 

clear all;
clc;

fid = fopen('plhs.dat','r');

fidp = fopen('sortie_pfem_py.txt', 'r');
fidm = fopen('sortie_mfem_py.txt', 'r');



A = dlmread('plhs.dat', ' ');

D1 = dlmread('sortie_pfem_py.txt', ' ');
D2 = dlmread('sortie_mfem_py.txt', ' ');


plot(A(1:100,3), D1(:,1),'+', A(1:100,3), D2(:,1),'*')


fclose(fid);
fclose(fidp);
fclose(fidm);