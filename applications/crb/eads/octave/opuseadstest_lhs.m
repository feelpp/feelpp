# -*- octave -*- 

clear all;
clc;

fid = fopen('mplhs.dat','r');

fidp = fopen('sortie_pfem_oct.txt', 'w');
fidm = fopen('sortie_mfem_oct.txt', 'w');

#############################
#### Default input values####
#############################

#  inP(1) = 10;   # kIC : thermal conductivity (default: 2)
#  inP(2) = 7e-3; # D : fluid flow rate (default: 5e-3)
#  inP(3) = 1e6;  # Q : heat flux (default: 1e6)
#  inP(4) = 100;  # r : conductance (default: 100)
#  inP(5) = 4e-3; # ea : length air flow channel (default: 4e-3)
#  inP(6) = 1; # meshsize times 1e-3 (default: 1)
#  inP(7) = 2; # integer : polynomial degree for the temperature in {1,2,3,4}

n=10;

A = dlmread('mplhs.dat', ' ');

for i=1:n
  inP(1) = A(i,1);
  inP(2) = 7.0e-3; #A(i,2);
  inP(3) = 1.0e+6;
  inP(4) = 100;  #A(i,3);
  inP(5) = 4.0e-3; #float(A[i,3])
  inP(6) = 1; # meshsize times 1e-3 (default: 1)
  inP(7) = 2; # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)

YP(i,:)= opuseadspfem( inP );

  inP(6) = 1; # meshsize times 1e-3 (default: 1)
  inP(7) = 2; # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
YM(i,:) = opuseadsmfem( inP );

end

dlmwrite('sortie_pfem_oct.txt', YP, ' ');
dlmwrite('sortie_mfem_oct.txt', YM, ' ');

fclose(fid);
fclose(fidp);
fclose(fidm);

#axis([0 120 320 330]);
#plot(D,y(1:2:10),'--og', D,y(2:2:10),'--rs');
