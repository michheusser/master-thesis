close all
clear all
clc

%% Numerical (Applied B)
load('SimulationB.mat')

for i = 1 : 14
   
    M = SimulationB.(['B' num2str(i)]).M;
    H = SimulationB.(['B' num2str(i)]).H;
    H_app = SimulationB.(['B' num2str(i)]).H_app;
    n_helix = SimulationB.(['B' num2str(i)]).n_helix;
    
    N_simulationB = (H_app - H)*inv(M);
    DemagFactors.(['B' num2str(i)]).N_simulationB = N_simulationB;

end


save('DemagFactors.mat','DemagFactors')
