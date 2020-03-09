%%
close all
clear all
clc
%%
mphstart()

%% REMANENT M
model = mphload('Helix1-10_(Remanent M).mph')
%mphnavigator
%info = mphsolinfo(model)
%%
for i = 1 : 10
    [H_x, H_y, H_z, M_x, M_y, M_z, M_app, n_helix_vec] = mphmean(model,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'M_app', 'n_helix'},'volume','selection',2,'dataset','dset2','outersolnum',i);
    %
    H = [H_x' H_y' H_z']';
    M = [M_x' M_y' M_z']';
    H_app = zeros(3,3);
    
    SimulationM.(['H' num2str(i)]).H = H;
    SimulationM.(['H' num2str(i)]).M = M;
    SimulationM.(['H' num2str(i)]).H_app = H_app;
    SimulationM.(['H' num2str(i)]).n_helix = n_helix_vec(1);
    SimulationM.(['H' num2str(i)]).Geometry = 'Helix';
end

save('SimulationM.mat','SimulationM')

%% APPLIED B FIELD
model = mphload('Helix1-10_(Applied B).mph')

mu_0 = 1.25663706e-6;
mu_r = 1;
mu = mu_r*mu_0;

for i = 1 : 10
    
    [H_x, H_y, H_z, M_x, M_y, M_z, B_app, n_helix_vec] = mphmean(model,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'B_app', 'n_helix'},'volume','selection',2,'dataset','dset2','outersolnum',i);
    
    H = [H_x' H_y' H_z']';
    M = [M_x' M_y' M_z']';
    H_app = B_app(1)/mu*eye(3);
    
    
    SimulationB.(['H' num2str(i)]).H = H;
    SimulationB.(['H' num2str(i)]).M = M;
    SimulationB.(['H' num2str(i)]).H_app = H_app;
    SimulationB.(['H' num2str(i)]).n_helix = n_helix_vec(1);
    SimulationB.(['H' num2str(i)]).Geometry = 'Helix';
    
end

save('SimulationB.mat','SimulationB')
