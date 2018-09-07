%# -*-matlab-*-

% Script for the SA using plhs sampling... 
% Second order sobol sensitivity indices

clear all;
clc;


fid   = fopen('plhs.dat', 'r');
fidp  = fopen('output_plhs_pfem_all_py.txt', 'r');

fidKicD = fopen('output_plhs_pfem_Kic_D_py.txt', 'r');
fidKicr  = fopen('output_plhs_pfem_Kic_r_py.txt', 'r');
fidKicea = fopen('output_plhs_pfem_Kic_ea_py.txt', 'r');
fidDr = fopen('output_plhs_pfem_D_r_py.txt', 'r');
fidDea = fopen('output_plhs_pfem_D_ea_py.txt', 'r');
fidrea = fopen('output_plhs_pfem_r_ea_py.txt', 'r');

% Sampling size
N=1000;

% P = [Kic D Q r ea]
% Kicmin = 2.0e-2; Kicmax = 1.5e+2;
% Dmin = 5.0e-4; Dmax = 1.0e-2;
% Q = 1.0e+6;
% Rmin = 0.1; Rmax = 1.0e+2;
% Eamin = 2.5e-3; Eamax = 5.0e-2;

% Read input matrix
P = dlmread('plhs.dat', ' ');

% Outputs related to input data: 'plhs.dat' all parameters are varied
S_pfem = dlmread('output_plhs_pfem_all_py.txt', ' ');

% Outputs related to input data: 'plhs.dat'
KicD = dlmread('output_plhs_pfem_Kic_D_py.txt', ' ');   % Kic and D are constants
Kicr  = dlmread('output_plhs_pfem_Kic_r_py.txt', ' ');  % Kic and r are constants
Kicea = dlmread('output_plhs_pfem_Kic_ea_py.txt', ' '); % Kic and ea are constants
Dr = dlmread('output_plhs_pfem_D_r_py.txt', ' ');       % D and r are constants
Dea  = dlmread('output_plhs_pfem_D_ea_py.txt', ' ');    % D and ea are constants
rea = dlmread('output_plhs_pfem_r_ea_py.txt', ' ');     % r and ea are constants


%%% Second Sobol Order Sensitivity Indices
for i=1:2
    Sob_ind_KicD(i) = var(KicD(:,i))/var(S_pfem(:,i));
    Sob_ind_Kicr(i) = var(Kicr(:,i))/var(S_pfem(:,i));
    Sob_ind_Kicea(i) = var(Kicea(:,i))/var(S_pfem(:,i));
    Sob_ind_Dr(i) = var(Dr(:,i))/var(S_pfem(:,i));
    Sob_ind_Dea(i) = var(Dea(:,i))/var(S_pfem(:,i));
    Sob_ind_rea(i) = var(rea(:,i))/var(S_pfem(:,i));
end

indice_sobol=[Sob_ind_KicD; Sob_ind_Kicr; Sob_ind_Kicea; Sob_ind_Dr; Sob_ind_Dea; Sob_ind_rea];


fclose(fid);
fclose(fidp);
fclose(fidKicD);
fclose(fidKicr);
fclose(fidKicea);
fclose(fidDr);
fclose(fidDea);
fclose(fidrea);

