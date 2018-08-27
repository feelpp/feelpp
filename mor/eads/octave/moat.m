%# -*-matlab-*-

% Script for the SA using the Morris One At a Time method...

clear all;
clc;


fid   = fopen('plhs.dat', 'r');

% Outputs related to input data: 'plhs.dat' all parameters are varied
fidp  = fopen('output_plhs_pfem_all_py.txt', 'r');

% Outputs related to input data: 'plhs.dat'
fidKic = fopen('output_plhs_pfem_Kic_py.txt', 'r');
fidD  = fopen('output_plhs_pfem_D_py.txt', 'r');
fidQ  = fopen('output_plhs_pfem_Q_py.txt', 'r');
fidr  = fopen('output_plhs_pfem_r_py.txt', 'r');
fidea = fopen('output_plhs_pfem_ea_py.txt', 'r');


% P = [Kic D Q r ea]
% Kicmin = 2.0e-2; Kicmax = 1.5e+2;
% Dmin = 5.0e-4; Dmax = 1.0e-2;
% Qmin = 0; Qmax = 1.0e+6;
% rmin = 0.1; rmax = 1.0e+2;
% eamin = 2.5e-3; eamax = 5.0e-2;

% Read input matrix
P = dlmread('plhs.dat', ' '); 

% Samlping size
N = 1000;

% Number of parameters
p=5;

S_pfem = dlmread('output_plhs_pfem_all_py.txt', ' ');

% Output for nominated input values
S = [326.570, 310.519];  

% Outputs related to input data: 'plhs.dat'
Kic = dlmread('output_plhs_pfem_Kic_py.txt', ' ');
D  = dlmread('output_plhs_pfem_D_py.txt', ' ');
Q  = dlmread('output_plhs_pfem_Q_py.txt', ' ');
r  = dlmread('output_plhs_pfem_r_py.txt', ' ');
ea = dlmread('output_plhs_pfem_ea_py.txt', ' ');

%%% Mean and standard deviation of input parameters
moy_P = mean(P);
std_P = std(P);

%%% Mean and standard deviation of outputs (S_1, S_2)
moy_pfem = mean(S_pfem);
std_pfem = std(S_pfem);

Reg_pfem(1) = regstats(S_pfem(:,1), P, 'linear');
Reg_pfem(2) = regstats(S_pfem(:,2), P, 'linear');

%%% Ridge regression coefficients (X^TX + kI) %%% 
 Ridge_pfem(:,1) = ridge(S_pfem(:,1), P, 1);
 Ridge_pfem(:,2) = ridge(S_pfem(:,2), P, 1);

%%% Correlation analysis %%%
for i=1:p
    Corre_S1_pfem(i) = corr(P(:,i),S_pfem(:,1));
    Corre_S2_pfem(i) = corr(P(:,i),S_pfem(:,2));
end

%%% S_1 Output
figure(1)
subplot(3,2,1);
scatter( P(:,1), Kic(:,1), 'r' )
Title=('"A"'); % ('Scatter plot pfem (S_1)');
title(Title);
xlabel('Conductivity K_{IC}');
ylabel('Temperature in K');

subplot(3,2,2);
scatter( P(:,2), D(:,1), 'r' )
Title=('"B"'); % ('Scatter plot pfem (S_1)');
title(Title);
xlabel('Inflow rate D');
ylabel('Temperature in K');

subplot(3,2,3);
scatter( P(:,3), Q(:,1), 'r' )
Title=('"C"'); % ('Scatter plot pfem (S_1)');
title(Title);
xlabel('Heat source Q');
ylabel('Temperature in K');

subplot(3,2,4);
scatter( P(:,4), r(:,1), 'r' )
Title=('"D"'); % ('Scatter plot pfem (S_1)');
title(Title);
xlabel('Thermal conductance r');
ylabel('Temperature in K');

subplot(3,2,5);
scatter( P(:,5), ea(:,1), 'r' )
Title=('"E"'); % ('Scatter plot pfem (S_1)');
title(Title);
xlabel('Thikness e_a');
ylabel('Temperature in K');


%%% S_2 Output
figure(2)
subplot(3,2,1);
scatter( P(:,1), Kic(:,2), 'r' )
Title=('"A"'); % ('Scatter plot pfem (S_2)');
title(Title);
xlabel('Conductivity K_{IC}');
ylabel('Temperature in K');

subplot(3,2,2);
scatter( P(:,2), D(:,2), 'r' )
Title=('"B"'); % ('Scatter plot pfem (S_2)');
title(Title);
xlabel('Inflow rate D');
ylabel('Temperature in K');

subplot(3,2,3);
scatter( P(:,3), Q(:,2), 'r' )
Title=('"C"'); % ('Scatter plot pfem (S_2)');
title(Title);
xlabel('Heat source Q');
ylabel('Temperature in K');

subplot(3,2,4);
scatter( P(:,4), r(:,2), 'r' )
Title=('"D"'); % ('Scatter plot pfem (S_2)');
title(Title);
xlabel('Thermal conductance r');
ylabel('Temperature in K');

subplot(3,2,5);
scatter( P(:,5), ea(:,2), 'r' )
Title=('"E"'); % ('Scatter plot pfem (S_2)');
title(Title);
xlabel('Thikness e_a');
ylabel('Temperature in K');


%%%% Sensitivity indices for the Morris OAT method

delta_Kic = P(:,1)-2.0e+0;
delta_S(:,1) = Kic(:,1) - S(1);
delta_S(:,2) = Kic(:,2) - S(2);
for i=1:N
    d(i,1) = delta_S(i,1)/delta_Kic(i);
    d(i,2) = delta_S(i,2)/delta_Kic(i);
end
indice_Kic = sum(abs(d))/N;

delta_D = P(:,2)-7.0e-3;
delta_S(:,1) = D(:,1) - S(1);
delta_S(:,2) = D(:,2) - S(2);
for i=1:1000
    d(i,1) = delta_S(i,1)/delta_D(i);
    d(i,2) = delta_S(i,2)/delta_D(i);
end
indice_D = sum(abs(d))/N;

delta_Q = P(:,3)-1.0e+6;
delta_S(:,1) = Q(:,1) - S(1);
delta_S(:,2) = Q(:,2) - S(2);
for i=1:N
    d(i,1) = delta_S(i,1)/delta_Q(i);
    d(i,2) = delta_S(i,2)/delta_Q(i);
end
indice_Q = sum(abs(d))/N;

delta_r = P(:,4)-1.0e+2;
delta_S(:,1) = r(:,1) - S(1);
delta_S(:,2) = r(:,2) - S(2);
for i=1:N
    d(i,1) = delta_S(i,1)/delta_r(i);
    d(i,2) = delta_S(i,2)/delta_r(i);
end
indice_r = sum(abs(d))/N;

delta_ea = P(:,5)-4.0e-3;
delta_S(:,1) = ea(:,1) - S(1);
delta_S(:,2) = ea(:,2) - S(2);
for i=1:N
    d(i,1) = delta_S(i,1)/delta_ea(i);
    d(i,2) = delta_S(i,2)/delta_ea(i);
end
indice_ea = sum(abs(d))/N;


indice_OAT = [indice_Kic; indice_D; indice_Q; indice_r; indice_ea];

fclose(fid);
fclose(fidp);
fclose(fidKic);
fclose(fidD);
fclose(fidQ);
fclose(fidr);
fclose(fidea);

