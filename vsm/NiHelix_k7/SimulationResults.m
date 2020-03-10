%%%%
close all
clear all
clc
%%
%mphstart()

%% High Field
model_High = mphload('Helix_VSM_k7_Nickel_High.mph')
%mphnavigator
%info = mphsolinfo(model)
%%
    [H_x, H_y, H_z, M_x, M_y, M_z, M_app] = mphmean(model_High,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'M_app'},'volume','selection',2,'dataset','dset1');
    H = [H_x' H_y' H_z']';
    M = [M_x' M_y' M_z']';
    H_app = zeros(3,3);
    N_simulation_High = (H_app - H)/M;
    
    SimulationResultsHigh.H = H;
    SimulationResultsHigh.M = M;
    SimulationResultsHigh.H_app = H_app;

save('SimulationResultsHigh.mat','SimulationResultsHigh')
save('N_simulation_High', 'N_simulation_High')

%% LowField
model_Low = mphload('Helix_VSM_k7_Nickel_Low.mph')
%%
mu_0 = 1.25663706e-6;
mu_r = 1;
mu = mu_r*mu_0;
    
    [H_x, H_y, H_z, M_x, M_y, M_z, B_app] = mphmean(model_Low,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'B_app'},'volume','selection',2,'dataset','dset1');
    
    H = [H_x' H_y' H_z']';
    M = [M_x' M_y' M_z']';
    H_app = B_app(1)/mu*eye(3);
    
    N_simulation_Low = (H_app - H)/M;
    
    SimulationResultsLow.H = H;
    SimulationResultsLow.M = M;
    SimulationResultsLow.H_app = H_app;

    
save('SimulationResultsLow.mat','SimulationResultsLow')
save('N_simulation_Low', 'N_simulation_Low')

%%
k_helix = 7
N_simulation_High
N_simulation_Low