%%
close all
clear all
clc
%%
mphstart()

%% REMANENT M
model = mphload('Helix1-10_ThinCoated(Remanent M).mph')
mphnavigator
%info = mphsolinfo(model)
%%
Helix_vec = [1:0.1:3 4:10]';
for i = 1 : length(Helix_vec)
    Helix_vec(i)
    [H_x, H_y, H_z, M_x, M_y, M_z, M_app, n_helix_vec] = mphmean(model,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'M_app', 'n_helix'},'volume','selection',2,'dataset','dset2','outersolnum',i);
    %
    H = [H_x' H_y' H_z']';
    M = [M_x' M_y' M_z']';
    H_app = zeros(3,3);
    
    SimulationM.(['S' num2str(i)]).H = H;
    SimulationM.(['S' num2str(i)]).M = M;
    SimulationM.(['S' num2str(i)]).H_app = H_app;
    SimulationM.(['S' num2str(i)]).n_helix = n_helix_vec(1);
    SimulationM.(['S' num2str(i)]).Geometry = 'Helix';
end

save('SimulationM.mat','SimulationM')

%% Numerical (Remanent M)


for i = 1 : length(Helix_vec)
    
    M = SimulationM.(['S' num2str(i)]).M;
    H = SimulationM.(['S' num2str(i)]).H;
    H_app = SimulationM.(['S' num2str(i)]).H_app;
    n_helix = SimulationM.(['S' num2str(i)]).n_helix;
    N_simulationM = (H_app - H)*inv(M);
    
    DemagFactors.(['S' num2str(i)]).N_simulationM = N_simulationM;
    DemagFactors.(['S' num2str(i)]).n_helix = Helix_vec(i);
end

save('DemagFactors.mat','DemagFactors')

