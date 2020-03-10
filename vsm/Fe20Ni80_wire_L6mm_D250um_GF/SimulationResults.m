%%%%
close all
clear all
clc
%%
%mphstart()

%% High Field
model_High = mphload('Helix_VSM_Wire_NiFe_HIGH.mph')
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

%% LowField 5.5mm
model_Low_55 = mphload('Helix_VSM_Wire_NiFe_Low_55.mph')
%%
mu_0 = 1.25663706e-6;
mu_r = 1;
mu = mu_r*mu_0;
    
    [H_x, H_y, H_z, M_x, M_y, M_z, B_app] = mphmean(model_Low_55,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'B_app'},'volume','selection',2,'dataset','dset1');
    
    H = [H_x' H_y' H_z']';
    M = [M_x' M_y' M_z']';
    H_app = B_app(1)/mu*eye(3);
    
    N_simulation_Low55 = (H_app - H)/M;
    
    SimulationResultsLow55.H = H;
    SimulationResultsLow55.M = M;
    SimulationResultsLow55.H_app = H_app;

    
save('SimulationResultsLow55.mat','SimulationResultsLow55')
save('N_simulation_Low55', 'N_simulation_Low55')


%% LowField 4.5mm
model_Low_45 = mphload('Helix_VSM_Wire_NiFe_Low_45.mph')
%%
mu_0 = 1.25663706e-6;
mu_r = 1;
mu = mu_r*mu_0;
    
    [H_x, H_y, H_z, M_x, M_y, M_z, B_app] = mphmean(model_Low_45,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'B_app'},'volume','selection',2,'dataset','dset1');
    
    H = [H_x' H_y' H_z']';
    M = [M_x' M_y' M_z']';
    H_app = B_app(1)/mu*eye(3);
    
    N_simulation_Low45 = (H_app - H)/M;
    
    SimulationResultsLow45.H = H;
    SimulationResultsLow45.M = M;
    SimulationResultsLow45.H_app = H_app;

    
save('SimulationResultsLow45.mat','SimulationResultsLow45')
save('N_simulation_Low45', 'N_simulation_Low45')

%% LowField 6.5mm
model_Low_65 = mphload('Helix_VSM_Wire_NiFe_Low_65.mph')
%%
mu_0 = 1.25663706e-6;
mu_r = 1;
mu = mu_r*mu_0;
    
    [H_x, H_y, H_z, M_x, M_y, M_z, B_app] = mphmean(model_Low_65,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'B_app'},'volume','selection',2,'dataset','dset1');
    
    H = [H_x' H_y' H_z']';
    M = [M_x' M_y' M_z']';
    H_app = B_app(1)/mu*eye(3);
    
    N_simulation_Low65 = (H_app - H)/M;
    
    SimulationResultsLow65.H = H;
    SimulationResultsLow65.M = M;
    SimulationResultsLow65.H_app = H_app;

    
save('SimulationResultsLow65.mat','SimulationResultsLow65')
save('N_simulation_Low65', 'N_simulation_Low65')


