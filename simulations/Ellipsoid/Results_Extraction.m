%%
close all
clear all
clc
%%
mphstart()

%% REMANENT M
model = mphload('Ellipsoid Saturated (Remanent M).mph')
%mphnavigator
%info = mphsolinfo(model)
[H_x, H_y, H_z, M_x, M_y, M_z, M_app, n_helix_vec] = mphmean(model,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'M_app', 'n_helix'},'volume','selection',2);

H = [H_x' H_y' H_z']';
M = [M_x' M_y' M_z']';
H_app = zeros(3,3);


SimulationM.('S1').H = H;
SimulationM.('S1').M = M;
SimulationM.('S1').H_app = H_app;
SimulationM.('S1').n_helix = n_helix_vec(1);
SimulationM.('S1').Geometry = 'Ellipsoid';

save('SimulationM.mat','SimulationM')

%% APPLIED B FIELD
model = mphload('Ellipsoid Saturated (Applied B).mph')

mu_0 = 1.25663706e-6;
mu_r = 1;
mu = mu_r*mu_0;

[H_x, H_y, H_z, M_x, M_y, M_z, B_app, n_helix_vec] = mphmean(model,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'B_app', 'n_helix'},'volume','selection',2);

H = [H_x' H_y' H_z']';
M = [M_x' M_y' M_z']';
H_app = B_app(1)/mu*eye(3);


SimulationB.('S1').H = H;
SimulationB.('S1').M = M;
SimulationB.('S1').H_app = H_app;
SimulationB.('S1').n_helix = n_helix_vec(1);
SimulationB.('S1').Geometry = 'Ellipsoid';


save('SimulationB.mat','SimulationB')
