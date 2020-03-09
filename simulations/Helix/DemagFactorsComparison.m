close all
clear all
clc
%%
Data_Full = load('Helix1-10(Detailed)_Saturated/DemagFactors.mat');
Data_ThinCoated = load('Helix1-10_Saturated_Thin_Coated/DemagFactors.mat');
Data_HalfCoated = load('Helix1-10_Saturated_Half_Coated/DemagFactors.mat');
Data_Full_Line_Analytical = load('DemagMatrix_Helix_Line_Average_Analytical.mat');
Data_Full_Vol_Analytical = load('DemagMatrix_Helix_Vol_Average_Analytical.mat');
Data_Full_Low_Vol_Simulation = load('Helix_Low_Fields/DemagFactors.mat')
Data_Thin_Low_Vol_Simulation = load('Helix_Thin_Low_Fields/DemagFactors.mat')
Data_Half_Low_Vol_Simulation = load('Helix_Half_Low_Fields/DemagFactors.mat')

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

N_11_Full_vec = zeros(length(Helix_full_vec),1);
N_22_Full_vec = zeros(length(Helix_full_vec),1);
N_33_Full_vec = zeros(length(Helix_full_vec),1);

N_11_ThinCoated_vec = zeros(length(Helix_thin_vec),1);
N_22_ThinCoated_vec = zeros(length(Helix_thin_vec),1);
N_33_ThinCoated_vec = zeros(length(Helix_thin_vec),1);

N_11_HalfCoated_vec = zeros(length(Helix_half_vec),1);
N_22_HalfCoated_vec = zeros(length(Helix_half_vec),1);
N_33_HalfCoated_vec = zeros(length(Helix_half_vec),1);

N_11_Full_Line_An_vec = zeros(length(Helix_full_line_an_vec),1);
N_22_Full_Line_An_vec = zeros(length(Helix_full_line_an_vec),1);
N_33_Full_Line_An_vec = zeros(length(Helix_full_line_an_vec),1);

N_11_Full_Vol_An_vec = zeros(length(Helix_full_vol_an_vec),1);
N_22_Full_Vol_An_vec = zeros(length(Helix_full_vol_an_vec),1);
N_33_Full_Vol_An_vec = zeros(length(Helix_full_vol_an_vec),1);

N_11_Low_Full_Vol_Sim_vec = zeros(length(Helix_low_full_vol_sim_vec),1);
N_22_Low_Full_Vol_Sim_vec = zeros(length(Helix_low_full_vol_sim_vec),1);
N_33_Low_Full_Vol_Sim_vec = zeros(length(Helix_low_full_vol_sim_vec),1);

N_11_Low_Thin_Vol_Sim_vec = zeros(length(Helix_low_thin_vol_sim_vec),1);
N_22_Low_Thin_Vol_Sim_vec = zeros(length(Helix_low_thin_vol_sim_vec),1);
N_33_Low_Thin_Vol_Sim_vec = zeros(length(Helix_low_thin_vol_sim_vec),1);

N_11_Low_Half_Vol_Sim_vec = zeros(length(Helix_low_half_vol_sim_vec),1);
N_22_Low_Half_Vol_Sim_vec = zeros(length(Helix_low_half_vol_sim_vec),1);
N_33_Low_Half_Vol_Sim_vec = zeros(length(Helix_low_half_vol_sim_vec),1);

for i = 1 : length(Helix_thin_vec)
    
    N_11_ThinCoated_vec(i) = DemagFactors_ThinCoated.(['S' num2str(i)]).N_simulationM(1,1);
    N_22_ThinCoated_vec(i) = DemagFactors_ThinCoated.(['S' num2str(i)]).N_simulationM(2,2);
    N_33_ThinCoated_vec(i) = DemagFactors_ThinCoated.(['S' num2str(i)]).N_simulationM(3,3);
   
    
end


for i = 1: length(Helix_half_vec)
    
    N_11_HalfCoated_vec(i) = DemagFactors_HalfCoated.(['S' num2str(i)]).N_simulationM(1,1);
    N_22_HalfCoated_vec(i) = DemagFactors_HalfCoated.(['S' num2str(i)]).N_simulationM(2,2);
    N_33_HalfCoated_vec(i) = DemagFactors_HalfCoated.(['S' num2str(i)]).N_simulationM(3,3);

end



for i = 1 : length(Helix_full_vec)
    
    N_11_Full_vec(i) = DemagFactors_Full.(['S' num2str(i)]).N_simulationM(1,1);
    N_22_Full_vec(i) = DemagFactors_Full.(['S' num2str(i)]).N_simulationM(2,2);
    N_33_Full_vec(i) = DemagFactors_Full.(['S' num2str(i)]).N_simulationM(3,3);
    
end

for i = 1 : length(Helix_full_line_an_vec)
    
    N_11_Full_Line_An_vec(i) = DemagFactors_Full_Line_Analytical.(['S' num2str(i)]).N_matrix_global_line_av(1,1);
    N_22_Full_Line_An_vec(i) = DemagFactors_Full_Line_Analytical.(['S' num2str(i)]).N_matrix_global_line_av(2,2);
    N_33_Full_Line_An_vec(i) = DemagFactors_Full_Line_Analytical.(['S' num2str(i)]).N_matrix_global_line_av(3,3);
   
end
   
for i = 1 : length(Helix_full_line_an_vec)
    
    N_11_Full_Vol_An_vec(i) = DemagFactors_Full_Vol_Analytical.(['S' num2str(i)]).N_matrix_global_vol_av(1,1);
    N_22_Full_Vol_An_vec(i) = DemagFactors_Full_Vol_Analytical.(['S' num2str(i)]).N_matrix_global_vol_av(2,2);
    N_33_Full_Vol_An_vec(i) = DemagFactors_Full_Vol_Analytical.(['S' num2str(i)]).N_matrix_global_vol_av(3,3);
   
end


for i = 1 : length(Helix_low_full_vol_sim_vec)
    
    N_11_Low_Full_Vol_Sim_vec(i) = DemagFactors_Full_Low_Vol_Simulation.(['B' num2str(i)]).N_simulationB(1,1);
    N_22_Low_Full_Vol_Sim_vec(i) = DemagFactors_Full_Low_Vol_Simulation.(['B' num2str(i)]).N_simulationB(2,2);
    N_33_Low_Full_Vol_Sim_vec(i) = DemagFactors_Full_Low_Vol_Simulation.(['B' num2str(i)]).N_simulationB(3,3);
   
end

for i = 1 : length(Helix_low_thin_vol_sim_vec)
    
    N_11_Low_Thin_Vol_Sim_vec(i) = DemagFactors_Thin_Low_Vol_Simulation.(['B' num2str(i)]).N_simulationB(1,1);
    N_22_Low_Thin_Vol_Sim_vec(i) = DemagFactors_Thin_Low_Vol_Simulation.(['B' num2str(i)]).N_simulationB(2,2);
    N_33_Low_Thin_Vol_Sim_vec(i) = DemagFactors_Thin_Low_Vol_Simulation.(['B' num2str(i)]).N_simulationB(3,3);
   
end

for i = 1 : length(Helix_low_half_vol_sim_vec)
    
    N_11_Low_Half_Vol_Sim_vec(i) = DemagFactors_Half_Low_Vol_Simulation.(['B' num2str(i)]).N_simulationB(1,1);
    N_22_Low_Half_Vol_Sim_vec(i) = DemagFactors_Half_Low_Vol_Simulation.(['B' num2str(i)]).N_simulationB(2,2);
    N_33_Low_Half_Vol_Sim_vec(i) = DemagFactors_Half_Low_Vol_Simulation.(['B' num2str(i)]).N_simulationB(3,3);
   
end

%%
subplot(3,1,1)
plot(Helix_full_vec,N_11_Full_vec,'*-',Helix_thin_vec,N_11_ThinCoated_vec,'o-',Helix_half_vec, N_11_HalfCoated_vec,'x-',Helix_full_line_an_vec,N_11_Full_Line_An_vec,'+-',Helix_full_vol_an_vec,N_11_Full_Vol_An_vec,'s-',Helix_low_full_vol_sim_vec,N_11_Low_Full_Vol_Sim_vec,'d-',Helix_low_thin_vol_sim_vec,N_11_Low_Thin_Vol_Sim_vec,'^-',Helix_low_half_vol_sim_vec,N_11_Low_Half_Vol_Sim_vec,'v-')
grid on
xlabel('Helix')
ylabel('N')
%ylim([0 1])
title({'','Comparison of demagnetization tensors for different helices in global coordinates','','N_{11}'})
legend('Full Magnetic: Volume Average (Simulation)', 'Thin Coated: Volume Average (Simulation)','Half Coated: Volume Average (Simulation)','Full Magnetic: Line Av (Quasi-Analytical)','Full Magnetic: Volumen Av (Quasi-Analytical)','Full Magnetic: Volumen Av, Low Fields (Simulation)', 'Thin Coated: Volumen Av, Low Fields (Simulation)' ,'Half Coated: Volume Av, Low Fields (Simulation)','Location','Best')

subplot(3,1,2)
plot(Helix_full_vec,N_22_Full_vec,'*-',Helix_thin_vec,N_22_ThinCoated_vec,'o-',Helix_half_vec,N_22_HalfCoated_vec,'x-',Helix_full_line_an_vec,N_22_Full_Line_An_vec,'+-',Helix_full_vol_an_vec,N_22_Full_Vol_An_vec,'s-',Helix_low_full_vol_sim_vec,N_22_Low_Full_Vol_Sim_vec,'d-',Helix_low_thin_vol_sim_vec,N_22_Low_Thin_Vol_Sim_vec,'^-',Helix_low_half_vol_sim_vec,N_22_Low_Half_Vol_Sim_vec,'v-')
grid on
xlabel('Helix')
ylabel('N')
%ylim([0 1])
title('N_{22}')
legend('Full Magnetic: Volume Average (Simulation)', 'Thin Coated: Volume Average (Simulation)','Half Coated: Volume Average (Simulation)','Full Magnetic: Line Av (Quasi-Analytical)','Full Magnetic: Volumen Av (Quasi-Analytical)','Full Magnetic: Volumen Av, Low Fields (Simulation)', 'Thin Coated: Volumen Av, Low Fields (Simulation)' ,'Half Coated: Volume Av, Low Fields (Simulation)','Location','Best')

subplot(3,1,3)
plot(Helix_full_vec,N_33_Full_vec,'*-',Helix_thin_vec,N_33_ThinCoated_vec,'o-',Helix_half_vec,N_33_HalfCoated_vec,'x-',Helix_full_line_an_vec,N_33_Full_Line_An_vec,'+-',Helix_full_vol_an_vec,N_33_Full_Vol_An_vec,'s-',Helix_low_full_vol_sim_vec,N_33_Low_Full_Vol_Sim_vec,'d-',Helix_low_thin_vol_sim_vec,N_33_Low_Thin_Vol_Sim_vec,'^-',Helix_low_half_vol_sim_vec,N_33_Low_Half_Vol_Sim_vec,'v-')
grid on
xlabel('Helix')
ylabel('N')
%ylim([0 1])
title('N_{33}')
legend('Full Magnetic: Volume Average (Simulation)', 'Thin Coated: Volume Average (Simulation)','Half Coated: Volume Average (Simulation)','Full Magnetic: Line Av (Quasi-Analytical)','Full Magnetic: Volumen Av (Quasi-Analytical)','Full Magnetic: Volumen Av, Low Fields (Simulation)', 'Thin Coated: Volumen Av, Low Fields (Simulation)' ,'Half Coated: Volume Av, Low Fields (Simulation)','Location','Best')
