%# -*-matlab-*-

% Script for the SA using parametrized latin hypercube sampling...

clear all;
clc;

fid   = fopen('plhs.dat', 'r');

fidp  = fopen('output_plhs_pfem_all_py.txt', 'r');

fidKic = fopen('output_plhs_pfem_Kic_py.txt', 'r');
fidD  = fopen('output_plhs_pfem_D_py.txt', 'r');
fidr  = fopen('output_plhs_pfem_r_py.txt', 'r');
fidea = fopen('output_plhs_pfem_ea_py.txt', 'r');

%%% Samlping size
N = 1000;

% Number of parameters
p=4;

% P = [Kic D Q r ea]
% Kicmin = 2.0e-2; Kicmax = 1.5e+2;
% Dmin = 5.0e-4; Dmax = 1.0e-2;
% Q = 1.0e+6;
% rmin = 0.1; rmax = 1.0e+2;
% eamin = 2.5e-3; eamax = 5.0e-2;

%%% Read input matrix
P = dlmread('plhs.dat', ' '); % or plhs.dat 

%%% Outputs related to input data: 'plhs.dat' all parameters are varied
S_pfem = dlmread('output_plhs_pfem_all_py.txt', ' ');

%%% Outputs related to input data: 'plhs.dat'
Kic = dlmread('output_plhs_pfem_Kic_py.txt', ' '); % Kic = 2
D  = dlmread('output_plhs_pfem_D_py.txt', ' ');    % D = 7.0e-3
r  = dlmread('output_plhs_pfem_r_py.txt', ' ');    % r =1.0e+2
ea = dlmread('output_plhs_pfem_ea_py.txt', ' ');   % ea = 4.0e-3

%%% Mean and standard deviation of input parameters
moy_P = mean(P);
std_P = std(P);

%%% Mean and standard deviation of outputs (S_1, S_2)
moy_pfem = mean(S_pfem);
std_pfem = std(S_pfem);

%%% Histograms of S_1 and S_2
figure(1)
subplot(1,2,1);
hist(S_pfem(:,1))
title('"A"'); % ('S_1 (PFEM)');
xlabel('Temperature in K');
ylabel('Frequency');

subplot(1,2,2);
hist(S_pfem(:,2))
title('"B"'); % ('S_2 (PFEM)');
xlabel('Temperature in K');
ylabel('Frequency');

saveas(figure(1), 'hist')

%%% Cumulative Distribution Function of Output
figure(2)
subplot(1,2,1);
cdfplot(S_pfem(:,1))
hold on
T = min(S_pfem(:,1)):1:max(S_pfem(:,1));
F = cdf('normal', T, moy_pfem(1), std_pfem(1));
plot(T,F,'m')
legend('Simulated','Theoretical','Location','NW')
xlabel('Temperature in K');
title('Empirical CDF vs Data for pfem (S_1)')
hold off
[h_S1_pfem,p_S1_pfem,k_S1_pfem] = kstest2(S_pfem(:,1), F); % Kolmogorov-Smironov Test

subplot(1,2,2);
cdfplot(S_pfem(:,2))
hold on
T = min(S_pfem(:,2)):1:max(S_pfem(:,2));
F = cdf('normal', T,moy_pfem(2), std_pfem(2));
plot(T,F,'m')
legend('Simulated','Theoretical','Location','NW')
xlabel('Temperature in K');
title('Empirical CDF vs Data for pfem (S_2)')
hold off
[h_S2_pfem,p_S2_pfem,k_S2_pfem] = kstest2(S_pfem(:,2),F); % Kolmogorov-Smironov Test

saveas(figure(2), 'cdf')

[B_pfem1 Stats_pfem1] = robustfit(P(:,:), S_pfem(:,1));
[B_pfem2 Stats_pfem2] = robustfit(P(:,:), S_pfem(:,2)); 


% Reg_pfem(1) = regstats(S_pfem(:,1), P, 'linear');

%%% Ridge regression coefficients (X^TX + kI) %%% 

 Ridge_pfem(:,1) = ridge(S_pfem(:,1), P, 1);
 Ridge_pfem(:,2) = ridge(S_pfem(:,2), P, 1);
%%% Correlation analysis %%%
for i=1:p
    Corre_S1_pfem(i) = corr(P(:,i),S_pfem(:,1));
    Corre_S2_pfem(i) = corr(P(:,i),S_pfem(:,2));
end

%%% First Sobol Order Sensitivity Indices
for i=1:2
    Sob_ind_K(i) = var(Kic(:,i))/var(S_pfem(:,i));
    Sob_ind_D(i) = var(D(:,i))/var(S_pfem(:,i));
    Sob_ind_r(i) = var(r(:,i))/var(S_pfem(:,i));
    Sob_ind_ea(i) = var(ea(:,i))/var(S_pfem(:,i));
end

%%% Scatter plot of S_1
figure(3)
subplot(2,2,1);
scatter( P(:,1), S_pfem(:,1), 'r' )
lsline
Title='"A"'; % 'Scatter plot pfem (S_1)';
title(Title);
xlabel('Conductivity K_{IC}');
ylabel('Temperature in K');
legend;
refline

subplot(2,2,2);
scatter( P(:,2), S_pfem(:,1), 'r' )
lsline
Title='"B"'; % 'Scatter plot pfem (S_1)';
title(Title);
xlabel('Inflow rate D');
ylabel('Temperature in K');

subplot(2,2,3);
scatter( P(:,3), S_pfem(:,1), 'r' )
lsline
Title='"C"'; % 'Scatter plot pfem (S_1)';
title(Title);
xlabel('Thermal conductance r');
ylabel('Temperature in K');

subplot(2,2,4);
scatter( P(:,4), S_pfem(:,1), 'r' )
lsline
Title='"D"'; % 'Scatter plot pfem (S_1)';
title(Title);
xlabel('Thikness e_{a}');
ylabel('Temperature in K');

saveas(figure(3), 'scatter of S1');

%%% Scatter plot of S_2
figure(4)
subplot(2,2,1);
scatter( P(:,1), S_pfem(:,2), 'r' )
lsline
Title='"A"'; % 'Scatter plot pfem (S_2)';
title(Title);
xlabel('Conductivity K_{IC}');
ylabel('Temperature in K');
legend;
refline

subplot(2,2,2);
scatter( P(:,2), S_pfem(:,2), 'r' )
lsline
Title='"B"'; % 'Scatter plot pfem (S_2)';
title(Title);
xlabel('Inflow rate D');
ylabel('Temperature in K');

subplot(2,2,3);
scatter( P(:,3), S_pfem(:,2), 'r' )
lsline
Title='"C"'; % 'Scatter plot pfem (S_2)';
title(Title);
xlabel('Thermal conductance r');
ylabel('Temperature in K');

subplot(2,2,4);
scatter( P(:,4), S_pfem(:,2), 'r' )
lsline
Title='"D"'; % 'Scatter plot pfem (S_2)';
title(Title);
xlabel('Thikness e_{a}');
ylabel('Temperature in K');

saveas(figure(4), 'scatter of S2');

%%% Box plot
figure(5)
boxplot([S_pfem(:,1) S_pfem(:,2)]) %['S_1' 'S_2']
title('Box plot (PFEM)');
xlabel('Outputs S_1 and S_2');
ylabel('Temperature in K');


figure(6)
subplot(2,2,1);
scatter( P(:,2), Kic(:,1), 'r' )
Title='"A"'; % ('Scatter plot pfem (S_1) K_{IC} = 2');
title(Title);
xlabel('Inflow rate D');
ylabel('Temperature in K');

subplot(2,2,2);
scatter( P(:,3), Kic(:,1), 'r' )
Title='"B"'; % ('Scatter plot pfem (S_1) K_{IC} = 2');
title(Title);
xlabel('Thermal conductance r');
ylabel('Temperature in K');

subplot(2,2,3);
scatter( P(:,4), Kic(:,1), 'r' )
Title='"C"'; % ('Scatter plot pfem (S_1) K_{IC} = 2');
title(Title);
xlabel('Thikness e_a');
ylabel('Temperature in K');


figure(7)
subplot(2,2,1);
scatter( P(:,1), D(:,1), 'r' )
Title='"A"'; % ('Scatter plot pfem (S_1) D = 7.0\times 10^{-3}');
title(Title);
xlabel('Conductivity K_{IC}');
ylabel('Temperature in K');

subplot(2,2,2);
scatter( P(:,3), D(:,1), 'r' )
Title='"B"'; % ('Scatter plot pfem (S_1) D = 7.0\times 10^{-3}');
title(Title);
xlabel('Thermal conductance r');
ylabel('Temperature in K');

subplot(2,2,3);
scatter( P(:,4), D(:,1), 'r' )
Title='"C"'; % ('Scatter plot pfem (S_1) D = 7.0\times 10^{-3}');
title(Title);
xlabel('Thickness e_a');
ylabel('Temperature in K');



figure(8)
subplot(2,2,1);
scatter( P(:,1), r(:,1), 'r' )
Title='"A"'; % ('Scatter plot pfem (S_1) r = 10^2');
title(Title);
xlabel('Conductivity K_{IC}');
ylabel('Temperature in K');

subplot(2,2,2);
scatter( P(:,2), r(:,1), 'r' )
Title='"B"'; % ('Scatter plot pfem (S_1) r = 10^2');
title(Title);
xlabel('Inflow rate D');
ylabel('Temperature in K');

subplot(2,2,3);
scatter( P(:,4), r(:,1), 'r' )
Title='"C"'; % ('Scatter plot pfem (S_1) r = 10^2');
title(Title);
xlabel('Thickness e_a');
ylabel('Temperature in K');


figure(9)
subplot(2,2,1);
scatter( P(:,1), ea(:,1), 'r' )
Title='"A"'; % ('Scatter plot pfem (S_1) e_a = 4.0\times 10^{-3}');
title(Title);
xlabel('Conductivity K_{IC}');
ylabel('Temperature in K');

subplot(2,2,2);
scatter( P(:,2), ea(:,1), 'r' )
Title='"B"'; % ('Scatter plot pfem (S_1) e_a = 4.0\times 10^{-3}');
title(Title);
xlabel('Inflow rate D');
ylabel('Temperature in K');

subplot(2,2,3);
scatter( P(:,3), ea(:,1), 'r' )
Title='"C"'; % ('Scatter plot pfem (S_1) e_a = 4.0\times 10^{-3}');
title(Title);
xlabel('Thermal conductance r');
ylabel('Temperature in K');


fclose(fid);
fclose(fidp);
fclose(fidKic);
fclose(fidD);
fclose(fidr);
fclose(fidea);

