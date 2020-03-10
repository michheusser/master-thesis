close all
clear all 
clc

%% Parameters
k_helix = 3

% M4 Screw
mu_0 = 1.25663706e-6;
Oe_to_Am = 1e3/(4*pi);
emu_to_Am2 = 1e-3;
M_s = 1.1/mu_0%0.87/mu_0;

d_helix = 250e-6;
D_helix = 4e-3;
h_helix = 0.7e-3;
L = k_helix*sqrt(h_helix^2+(D_helix*pi)^2);
A = (d_helix/2)^2*pi;
V = L*A;

%%

Data_ax = dlmread('FeNiHelix1-VIRaxial.txt', ',', 1, 0);
H_app_ax = Data_ax(:,8)*Oe_to_Am;
M_ax_ax = (Data_ax(:,10)/1000)*emu_to_Am2/V;
M_ax_rad = (Data_ax(:,11)/1000)*emu_to_Am2/V;


M_s_ax_ax = max(M_ax_ax);
M_s_ax_rad = max(M_ax_ax);

Data_rad = dlmread('FeNiHelix1-VIRradial.txt', ',', 1, 0);
H_app_rad = Data_rad(:,8)*Oe_to_Am;
M_rad_rad = (Data_rad(:,10)/1000)*emu_to_Am2/V;
M_rad_ax = -(Data_rad(:,11)/1000)*emu_to_Am2/V;

M_s_rad_rad = max(M_rad_rad);
M_s_rad_ax = max(M_rad_ax);

plot(H_app_ax*mu_0,M_ax_ax,'-*',H_app_ax*mu_0,M_ax_rad,'-*',H_app_rad*mu_0,M_rad_ax,'-*',H_app_rad*mu_0,M_rad_rad,'-*')
grid on
axis tight
xlabel('Applied Field B [T]')
ylabel('Magnetization M [A/m]')
legend('M_{ax,ax}','M_{ax,rad}','M_{rad,ax}','M_{rad,rad}')

Helix_Data_Raw.B_app_ax = Data_ax(:,8);
Helix_Data_Raw.B_app_rad = Data_rad(:,8);
Helix_Data_Raw.M_ax_ax = Data_ax(:,10);
Helix_Data_Raw.M_ax_rad = Data_ax(:,11);
Helix_Data_Raw.M_rad_rad = Data_rad(:,10);
Helix_Data_Raw.M_rad_ax = Data_rad(:,11);

save('Helix_Data_Raw.mat','Helix_Data_Raw')

Helix_Data.B_app_ax = H_app_ax*mu_0;
Helix_Data.B_app_rad = H_app_rad*mu_0;
Helix_Data.M_ax_ax = M_ax_ax;
Helix_Data.M_ax_rad = M_ax_rad;
Helix_Data.M_rad_rad = M_rad_rad;
Helix_Data.M_rad_ax = M_rad_ax;

save('Helix_Data.mat','Helix_Data')

Helix_Data_Normalized.B_app_ax = H_app_ax*mu_0;
Helix_Data_Normalized.B_app_rad = H_app_rad*mu_0;
Helix_Data_Normalized.M_ax_ax = M_ax_ax/M_s_ax_ax;
Helix_Data_Normalized.M_ax_rad = M_ax_rad/M_s_ax_rad;
Helix_Data_Normalized.M_rad_rad = M_rad_rad/M_s_rad_rad;
Helix_Data_Normalized.M_rad_ax = M_rad_ax/M_rad_ax;

save('Helix_Data_Normalized.mat','Helix_Data_Normalized')

%% Simplification of Hysteresis Model
[H_app_ax_simplified, M_ax_ax_simplified, m_scale_ax] = SimplifyMH(H_app_ax,M_ax_ax,M_s,'Saturates');
[~, M_ax_rad_simplified ] = SimplifyMH(H_app_ax,M_ax_rad,m_scale_ax,'Zero');
[H_app_rad_simplified, M_rad_rad_simplified, m_scale_rad] = SimplifyMH(H_app_rad,M_rad_rad,M_s,'Saturates');
[~, M_rad_ax_simplified ] = SimplifyMH(H_app_rad,M_rad_ax,m_scale_rad,'Zero');

%% Splitting of Simplified Model (Positive Part)
H_app_monotone_ax = 0.22/mu_0;
I_pos_ax = find(H_app_ax_simplified >=0 & H_app_ax_simplified < H_app_monotone_ax);
H_app_ax_pos = H_app_ax_simplified(I_pos_ax);
M_ax_ax_pos = M_ax_ax_simplified(I_pos_ax);
M_ax_rad_pos = M_ax_rad_simplified(I_pos_ax);

H_app_monotone_rad = 1/mu_0;
I_pos_rad = find(H_app_rad_simplified >=0 & H_app_ax_simplified < H_app_monotone_rad);
H_app_rad_pos = H_app_rad_simplified(I_pos_rad);
M_rad_ax_pos = M_rad_ax_simplified(I_pos_rad);
M_rad_rad_pos = M_rad_rad_simplified(I_pos_rad);


%% Variance of Fit Analysis (High Fields)
N_end = 1;
[var_ax_high,M_s_ax_ax, M_s_ax_rad] = var_est_high(H_app_ax_pos,M_ax_ax_pos,M_ax_rad_pos,N_end);
[var_rad_high,M_s_rad_ax, M_s_rad_rad] = var_est_high(H_app_rad_pos,M_rad_ax_pos,M_rad_rad_pos,N_end);

figure(2)
plot(H_app_ax_pos*mu_0,var_ax_high,'-*',H_app_rad_pos*mu_0,var_rad_high,'-*')
grid on
axis tight
legend('Axial Direction', 'Radial Direction')
title({['Variance of linear fit for variable data points'],['Helix: k=' num2str(k_helix)]})
ylabel('Variance s^2 [A^2/m^2]')
xlabel('B_{app} [T] (Points Used)')

%% Variance of Fit Analysis (Low Fields)
[m_ax, var_ax_low] = var_est_low(H_app_ax_pos,M_ax_ax_pos,M_ax_rad_pos);
[m_rad, var_rad_low] = var_est_low(H_app_rad_pos,M_rad_ax_pos,M_rad_rad_pos);

figure(1)
plot(H_app_ax_pos*mu_0,var_ax_low,'-*',H_app_rad_pos*mu_0,var_rad_low,'-*')
grid on
axis tight
legend('Axial Direction', 'Radial Direction')
title({['Variance of linear fit for variable data points'],['Helix: k=' num2str(k_helix)]})
ylabel('Variance s^2 [A^2/m^2]')
xlabel('B_{app} [T] (Points Used)')

B_ax = 0.01; %Read from the Variance analysis plot
i_axial = length(find(H_app_ax_pos*mu_0 <= B_ax));

B_rad = 0.01; %Read from the Variance analysis plot
i_rad = length(find(H_app_rad_pos*mu_0 <= B_rad));

%% N simulation High
load('N_simulation_High.mat')
chi_a_sim_High = inv(N_simulation_High);

%% N simulation Low
load('N_simulation_Low.mat')
chi_a_sim_Low = inv(N_simulation_Low);

%% N measurements High
r_sat_ax = 0.89; %Saturation Point
r_sat_rad = 0.8724;

M_ax_ax_sat = r_sat_ax*M_s;
H_app_sat_ax  = interp1(M_ax_ax_pos,H_app_ax_pos,M_ax_ax_sat,'spline');
chi_a_High_ax_ax = M_s/H_app_sat_ax;
chi_a_High_ax_rad = 0;

M_rad_rad_sat = r_sat_rad*M_s;
H_app_sat_rad  = interp1(M_rad_rad_pos,H_app_rad_pos,M_rad_rad_sat,'spline');
chi_a_High_rad_rad = M_s/H_app_sat_rad;
chi_a_High_rad_ax = 0;

chi_a_High = [...
    chi_a_High_rad_rad   0                   chi_a_High_ax_rad;...
    0                   chi_a_High_rad_rad   chi_a_High_ax_rad;...
    chi_a_High_rad_ax    chi_a_High_rad_ax    chi_a_High_ax_ax];

N_high = inv(chi_a_High);

%% N measurements Low
chi_a_Low_ax_ax = m_ax(i_axial,1);
chi_a_Low_ax_rad = m_ax(i_axial,2);
chi_a_Low_rad_ax = m_rad(i_rad,1);
chi_a_Low_rad_rad = m_rad(i_rad,2);

chi_a_Low = [...
    chi_a_Low_rad_rad   0                   chi_a_Low_ax_rad;...
    0                   chi_a_Low_rad_rad   chi_a_Low_ax_rad;...
    chi_a_Low_rad_ax    chi_a_Low_rad_ax    chi_a_Low_ax_ax];

N_low = inv(chi_a_Low)

%% Plot of Curves with Fits
subplot(2,1,1)
plot(H_app_ax_pos*mu_0,M_ax_ax_pos,'-*',...
    H_app_ax_pos*mu_0,M_ax_rad_pos,'-*',...
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_Low(3,3),'b--',...
H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_sim_Low(3,3),'b-',...    
H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_Low(2,3),'r--',...
        H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_sim_Low(2,3),'r-',...
H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_High(3,3),'k--',...
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_sim_High(3,3),'k-',...
    H_app_ax_pos*mu_0,ones(length(H_app_ax_pos),1)*M_s)
axis tight
grid on
ylim([1.1*min([M_ax_ax_pos; M_ax_rad_pos]) 1.1*max([M_ax_ax_pos; M_ax_rad_pos])])
legend('Axial Magnetization','Radial Magnetization','Axial Linearization (Measurement)','Axial Linearization (Simulation)','Radial Linearization (Measurement)','Radial Linearization (Simulation)','Saturation Slope (Measurement)','Saturation Slope (Simulation)','Total Saturation')
title({['','Analysis of VSM Measurements (k = '  num2str(k_helix) ')'],'','Simplified Magnetization Curves (Axial Direction)'})
xlabel('Applied Field B [T]')
ylabel('Magnetization M [A/m]')

subplot(2,1,2)
plot(H_app_rad_pos*mu_0,M_rad_rad_pos,'-*',...
    H_app_rad_pos*mu_0,M_rad_ax_pos,'-*',...
    H_app_rad_pos*mu_0,H_app_rad_pos*chi_a_Low(1,1),'b--',...
    H_app_rad_pos*mu_0,H_app_rad_pos*chi_a_sim_Low(1,1),'b-',...
    H_app_rad_pos*mu_0,H_app_rad_pos*chi_a_Low(3,2),'r--',...
    H_app_rad_pos*mu_0,H_app_rad_pos*chi_a_sim_Low(3,2),'r-',...
    H_app_rad_pos*mu_0,H_app_rad_pos*chi_a_High(1,1),'k--',...
    H_app_rad_pos*mu_0,H_app_rad_pos*chi_a_sim_High(1,1),'k-',...
    H_app_rad_pos*mu_0,ones(length(H_app_rad_pos),1)*M_s)

axis tight
grid on
ylim([1.1*min([M_ax_ax_pos; M_ax_rad_pos]) 1.1*max([M_ax_ax_pos; M_ax_rad_pos])])
legend('Radial Magnetization','Axial Magnetization','Radial Linearization (Measurement)','Radial Linearization (Simulation)','Axial Linearization (Measurement)','Axial Linearization (Simulation)','Saturation Slope (Measurement)','Saturation Slope (Simulation)','Total Saturation')
title({['','Analysis of VSM Measurements (k = '  num2str(k_helix) ')'],'','Simplified Magnetization Curves (Radial Direction)'})
xlabel('Applied Field B [T]')
ylabel('Magnetization M [A/m]')





%%
N_simulation_Low
N_low
N_simulation_High
N_high

chi_a_sim_Low
chi_a_Low
chi_a_sim_High
chi_a_High
