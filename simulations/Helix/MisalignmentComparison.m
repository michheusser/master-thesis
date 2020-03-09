close all
clear all
clc
%%
Data_Full = load('Helix1-10(Detailed)_Saturated/DemagFactors.mat');
Data_ThinCoated = load('Helix1-10_Saturated_Thin_Coated/DemagFactors.mat');
Data_HalfCoated = load('Helix1-10_Saturated_Half_Coated/DemagFactors_High.mat');
Data_Full_Line_Analytical = load('DemagMatrix_Helix_Line_Average_Analytical.mat');
Data_Full_Vol_Analytical = load('DemagMatrix_Helix_Vol_Average_Analytical.mat');
Data_Full_Low_Vol_Simulation = load('Helix_Low_Fields/DemagFactors.mat');
Data_Thin_Low_Vol_Simulation = load('Helix_Thin_Low_Fields/DemagFactors.mat');
Data_Half_Low_Vol_Simulation = load('Helix_Half_Low_Fields/DemagFactors_Low.mat');

DemagFactors_Full = Data_Full.DemagFactors;
DemagFactors_ThinCoated = Data_ThinCoated.DemagFactors;
DemagFactors_HalfCoated = Data_HalfCoated.DemagFactors;
DemagFactors_Full_Line_Analytical = Data_Full_Line_Analytical.DemagMatrix_Line_Average_Analytical;
DemagFactors_Full_Vol_Analytical = Data_Full_Vol_Analytical.DemagMatrix_Vol_Average_Analytical;
DemagFactors_Full_Low_Vol_Simulation = Data_Full_Low_Vol_Simulation.DemagFactors;
DemagFactors_Thin_Low_Vol_Simulation = Data_Thin_Low_Vol_Simulation.DemagFactors;
DemagFactors_Half_Low_Vol_Simulation = Data_Half_Low_Vol_Simulation.DemagFactors;

%%
Helix_thin_vec = [1:0.1:3 4:10]';
Helix_full_vec = [1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 4 5 6 7 8 9 10]';
Helix_half_vec = [1.5 1.65 1.75 1.9 2 2.2 2.3 2.35 2.45 2.55 2.6 2.7 2.8 3 3.2 3.3 4 5 6 7.2 8 9.1 10]';
Helix_full_line_an_vec = [1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 4 5 6 7 8 9 10]';
Helix_full_vol_an_vec = [1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 4 5 6 7 8 9 10]';
Helix_low_full_vol_sim_vec = [1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 4 5 6 7 8 9 10]';
Helix_low_thin_vol_sim_vec = [1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 4 5 6 7 8 9 10]';
Helix_low_half_vol_sim_vec = [1.5 1.65 1.75 1.9 2 2.2 2.3 2.35 2.45 2.55 2.6 2.7 2.8 3 3.2 3.3 4 5 6 7.2 8 9.1 10]';

theta_Full_vec = zeros(length(Helix_full_vec),1);
theta_ThinCoated_vec = zeros(length(Helix_thin_vec),1);
theta_HalfCoated_vec = zeros(length(Helix_half_vec),1);
theta_Full_Line_An_vec = zeros(length(Helix_full_line_an_vec),1);
theta_Full_Vol_An_vec = zeros(length(Helix_full_vol_an_vec),1);
theta_Low_Full_Vol_Sim_vec = zeros(length(Helix_low_full_vol_sim_vec),1);
theta_Low_Thin_Vol_Sim_vec = zeros(length(Helix_low_thin_vol_sim_vec),1);
theta_Low_Half_Vol_Sim_vec = zeros(length(Helix_low_half_vol_sim_vec),1);

for i = 1 : length(Helix_thin_vec)
    
    N = DemagFactors_ThinCoated.(['S' num2str(i)]).N_simulationM;
    [V,D] = eig(N);
    [~,I] = min(diag(D));
    v = V(:,I);
    theta_ThinCoated_vec(i) = acos(abs(v(3))/norm(v));
    
end

for i = 1: length(Helix_half_vec)
    
    N = DemagFactors_HalfCoated.(['S' num2str(i)]).N_simulationM;
    [V,D] = eig(N);
    [~,I] = min(diag(D));
    v = V(:,I);
    theta_HalfCoated_vec(i) = acos(abs(v(3))/norm(v));
end

for i = 1 : length(Helix_full_vec)
    
    N = DemagFactors_Full.(['S' num2str(i)]).N_simulationM;
    [V,D] = eig(N);
    [~,I] = min(diag(D));
    v = V(:,I);
    theta_Full_vec(i) = acos(abs(v(3))/norm(v));
    
end

for i = 1 : length(Helix_full_line_an_vec)
    
    N = DemagFactors_Full_Line_Analytical.(['S' num2str(i)]).N_matrix_global_line_av;
    [V,D] = eig(N);
    [~,I] = min(diag(D));
    v = V(:,I);
    theta_Full_Line_An_vec(i) = acos(abs(v(3))/norm(v));
    
end

for i = 1 : length(Helix_full_vol_an_vec)
    
    N = DemagFactors_Full_Vol_Analytical.(['S' num2str(i)]).N_matrix_global_vol_av;
    [V,D] = eig(N);
    [~,I] = min(diag(D));
    v = V(:,I);
    theta_Full_Vol_An_vec(i) = acos(abs(v(3))/norm(v));
    
end

for i = 1 : length(Helix_low_full_vol_sim_vec)
    
    N = DemagFactors_Full_Low_Vol_Simulation.(['B' num2str(i)]).N_simulationB;
    [V,D] = eig(N);
    [~,I] = min(diag(D));
    v = V(:,I);
    theta_Low_Full_Vol_Sim_vec(i) = acos(abs(v(3))/norm(v));
    
end

for i = 1 : length(Helix_low_thin_vol_sim_vec)
    
    N = DemagFactors_Thin_Low_Vol_Simulation.(['B' num2str(i)]).N_simulationB;
    [V,D] = eig(N);
    [~,I] = min(diag(D));
    v = V(:,I);
    theta_Low_Thin_Vol_Sim_vec(i) = acos(abs(v(3))/norm(v));
    
    
end

for i = 1 : length(Helix_low_half_vol_sim_vec)
    
    N = DemagFactors_Half_Low_Vol_Simulation.(['B' num2str(i)]).N_simulationB;
    [V,D] = eig(N);
    [~,I] = min(diag(D));
    v = V(:,I);
    theta_Low_Half_Vol_Sim_vec(i) = acos(abs(v(3))/norm(v));
    
end

%%

plot(Helix_full_vec,theta_Full_vec*360/(2*pi),'*-',Helix_thin_vec,theta_ThinCoated_vec*360/(2*pi),'o-',Helix_half_vec, theta_HalfCoated_vec*360/(2*pi),'x-',Helix_full_line_an_vec,theta_Full_Line_An_vec*360/(2*pi),'+-',Helix_full_vol_an_vec,theta_Full_Vol_An_vec*360/(2*pi),'s-',Helix_low_full_vol_sim_vec,theta_Low_Full_Vol_Sim_vec*360/(2*pi),'d-',Helix_low_thin_vol_sim_vec,theta_Low_Thin_Vol_Sim_vec*360/(2*pi),'^-',Helix_low_half_vol_sim_vec,theta_Low_Half_Vol_Sim_vec*360/(2*pi),'v-')
grid on
xlabel('Helix')
ylabel('\theta')
%ylim([0 1])
title({'','Comparison of misalignment angle towards z-Axis for different helices',''})
legend('Full Magnetic: Volume Average (Simulation)', 'Thin Coated: Volume Average (Simulation)','Half Coated: Volume Average (Simulation)','Full Magnetic: Line Av (Quasi-Analytical)','Full Magnetic: Volumen Av (Quasi-Analytical)','Full Magnetic: Volumen Av, Low Fields (Simulation)', 'Thin Coated: Volumen Av, Low Fields (Simulation)' ,'Half Coated: Volume Av, Low Fields (Simulation)','Location','Best')
%legend('Full Magnetic: Volume Average (Simulation)', 'Thin Coated: Volume Average (Simulation)','Half Coated: Volume Average (Simulation)','Full Magnetic: Line Av (Quasi-Analytical)', 'Thin Coated: Volumen Av, Low Fields (Simulation)' ,'Half Coated: Volume Av, Low Fields (Simulation)','Location','Best')
hold on

%Experimental Misalignment Angle

load('Experimental_Misalignment.mat')

MisAngle_Av = mean(MisAngle);
MisAngle_SD = sqrt(var(MisAngle,0,1));
errorbar(2:9,MisAngle_Av,MisAngle_SD)
hold off
