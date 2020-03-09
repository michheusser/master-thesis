%
close all
clear all
clc
%
mphstart()

%% APPLIED B FIELD
model = mphload('Helix3_LowFields (Applied Field).mph')
%%
mu_0 = 1.25663706e-6;
mu_r = 1;
mu = mu_r*mu_0;

for i = 1 : 14
    i
    [H_x, H_y, H_z, M_x, M_y, M_z, B_app, n_helix_vec] = mphmean(model,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'B_app', 'n_helix'},'volume','selection',2,'dataset','dset2','outersolnum',i);
    
    H = [H_x' H_y' H_z'];
    M = [M_x' M_y' M_z'];
    H_app = B_app(1)/mu*eye(3);
    
    
    SimulationB.(['B' num2str(i)]).H = H;
    SimulationB.(['B' num2str(i)]).M = M;
    SimulationB.(['B' num2str(i)]).H_app = H_app;
    SimulationB.(['B' num2str(i)]).B_app = B_app(1);
    SimulationB.(['B' num2str(i)]).n_helix = n_helix_vec(1);
    SimulationB.(['B' num2str(i)]).Geometry = 'Helix';
    
end

save('SimulationB.mat','SimulationB')

%% Numerical (Applied B)
%load('SimulationB.mat')

for i = 1 : 14
   
    M = SimulationB.(['B' num2str(i)]).M;
    H = SimulationB.(['B' num2str(i)]).H;
    H_app = SimulationB.(['B' num2str(i)]).H_app;
    n_helix = SimulationB.(['B' num2str(i)]).n_helix;
    
    N_simulationB = (H_app - H)*inv(M);
    DemagFactors.(['B' num2str(i)]).N_simulationB = N_simulationB;

end


save('DemagFactors.mat','DemagFactors')
