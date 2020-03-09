close all
clear all
clc

%Block
%% Numerical (Remanent M)
load('SimulationM.mat')

M = SimulationM.S1.M;
H = SimulationM.S1.H;
H_app = SimulationM.S1.H_app;
n_helix = SimulationM.S1.n_helix;

N_simulationM = (H_app - H)*inv(M);

%% Numerical (Applied B)
load('SimulationB.mat')

M = SimulationB.S1.M;
H = SimulationB.S1.H;
H_app = SimulationB.S1.H_app;
n_helix = SimulationB.S1.n_helix;

N_simulationB = (H_app - H)*inv(M);

%% Analytical
[h_helix, D_helix] = Helixinfo(n_helix);
k_helix = 3;

a = D_helix/2;
b = D_helix/2;
c = h_helix*k_helix/2;

N_analytical = Demagfactor_Block(a,b,c);

DemagFactors.('S1').N_simulationB = N_simulationB;
DemagFactors.('S1').N_simulationM = N_simulationM;
DemagFactors.('S1').N_analytical = N_analytical;
save('DemagFactors.mat','DemagFactors')
%%

N_simulationM
N_simulationB
N_analytical

Error1 = diag(abs(N_analytical - N_simulationM)./N_analytical)*100
Error2 = diag(abs(N_analytical - N_simulationB)./N_analytical)*100

% disp([num2str(Error(1)*100) '%'])
% disp([num2str(Error(2)*100) '%'])
% disp([num2str(Error(3)*100) '%'])