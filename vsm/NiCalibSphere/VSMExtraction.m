close all
clear all
clc

%% Parameters

mu_0 = 1.25663706e-6;
Oe_to_Am = 1e3/(4*pi);
emu_to_Am2 = 1e-3;

d_helix = 250e-6;
D_helix = 1.51e-3 + d_helix;
h_helix = 0.4e-3;
d_sphere = 2.383e-3; %m
V = (4/3)*pi*(d_sphere/2)^3;

rho_Ni = 8914e3; %g/m^3
m_s = 54.97;  %emu/g or mA * m^2/g
M_s = m_s*rho_Ni*emu_to_Am2;

%% Data Extraction
Data = dlmread('Clean/NisphereFineRes2-VIR-01.txt', ',', 1, 0);
H_app = Data(:,8)*Oe_to_Am;
M_data = (Data(:,10)/1000)*emu_to_Am2/V;
r_scaling = M_s/M_data(end);

M = M_data*r_scaling;


[H_app_simplified, M_simplified ] = SimplifyMH(H_app,M,M_s,'Saturates');


%%

I_pos = find(H_app_simplified >=0);
H_app_pos = H_app_simplified(I_pos);
M_pos = M_simplified(I_pos);

[m_vec, var_low] = var_est_low(H_app_pos,M_pos);

i = 50;
M_s = max(M_pos)
m_av = sum(m_vec(1:i))/i;
M_low = m_av*H_app_pos;
M_high = ones(length(H_app_pos),1)*M_s;

H_sat = (M_s/m_av);
M_sat = interp1(H_app_pos,M_pos,H_sat,'spline')

r_sat = M_sat/M_s

figure(1)
plot(H_app_pos*mu_0,M_pos,'-*',H_app_pos*mu_0,M_low,H_app_pos*mu_0,M_high)
grid on
axis tight
ylim([min(M_pos) max(M_pos)])
xlabel('Applied Field B [T]')
ylabel('Magnetization M [A/m]')
title('Magnetization Curve for a calibration sphere')

n = 20;
chi_a = M_pos./H_app_pos;

chi_m = 1./(- 1/3 + 1./chi_a);

plot(H_app_pos*mu_0,chi_m)
grid on
axis tight

save('chi_m_NiSph.mat','H_app_pos','chi_m')

%%


