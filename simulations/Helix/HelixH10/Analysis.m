close all
clear all
clc

%% Numerical (Remanent M)
load('SimulationM.mat')

for i = 1 : 10
    
    M = SimulationM.(['H' num2str(i)]).M;
    H = SimulationM.(['H' num2str(i)]).H;
    H_app = SimulationM.(['H' num2str(i)]).H_app;
    n_helix = SimulationM.(['H' num2str(i)]).n_helix;
    N_simulationM = (H_app - H)*inv(M);
    
    DemagFactors.(['H' num2str(i)]).N_simulationM = N_simulationM;
end

%% Numerical (Applied B)
load('SimulationB.mat')

for i = 1 : 10
   
    M = SimulationB.(['H' num2str(i)]).M;
    H = SimulationB.(['H' num2str(i)]).H;
    H_app = SimulationB.(['H' num2str(i)]).H_app;
    n_helix = SimulationB.(['H' num2str(i)]).n_helix;
    
    N_simulationB = (H_app - H)*inv(M);
    DemagFactors.(['H' num2str(i)]).N_simulationB = N_simulationB;

end


save('DemagFactors.mat','DemagFactors')
