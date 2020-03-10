close all
clear all
clc

%% Parameters


mu_0 = 1.25663706e-6;
Oe_to_Am = 1e3/(4*pi);
emu_to_Am2 = 1e-3;
rho_Ni = 8914e3; %g/m^3
m_s = 54.97;  %emu/g or mA * m^2/g
M_s = m_s*rho_Ni*emu_to_Am2;
chi_m = 24;


d_helix = 250e-6;
L = 5e-3;
A = (d_helix/2)^2*pi;
V = L*A;

%%

Data_ax = dlmread('Clean/NiWire_calib2_parallel-VIR-00.txt', ',', 1, 0);
H_app_ax = Data_ax(:,8)*Oe_to_Am;
M_ax_ax = (Data_ax(:,10)/1000)*emu_to_Am2/V;
M_ax_rad = (Data_ax(:,11)/1000)*emu_to_Am2/V;


M_s_ax_ax = max(M_ax_ax);
M_s_ax_rad = max(M_ax_ax);

Data_rad = dlmread('Clean/NiWire_calib2_perpendicular-VI-01.txt', ',', 1, 0);
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
I_pos_ax = find(H_app_ax_simplified > 0);
H_app_ax_pos = [0; H_app_ax_simplified(I_pos_ax)];
M_ax_ax_pos = [0; M_ax_ax_simplified(I_pos_ax)];
M_ax_rad_pos = [0; M_ax_rad_simplified(I_pos_ax)];

I_pos_rad = find(H_app_rad_simplified > 0);
H_app_rad_pos = [0; H_app_rad_simplified(I_pos_rad)];
M_rad_ax_pos = [0; M_rad_ax_simplified(I_pos_rad)];
M_rad_rad_pos = [0; M_rad_rad_simplified(I_pos_rad)];

%% Variance of Fit Analysis (High Fields)
N_end = 1;
[var_ax_high,M_s_ax_ax, M_s_ax_rad] = var_est_high(H_app_ax_pos,M_ax_ax_pos,M_ax_rad_pos,N_end);
[var_rad_high,M_s_rad_ax, M_s_rad_rad] = var_est_high(H_app_rad_pos,M_rad_ax_pos,M_rad_rad_pos,N_end);

figure(2)
plot(H_app_ax_pos*mu_0,var_ax_high,'-*',H_app_rad_pos*mu_0,var_rad_high,'-*')
grid on
axis tight
legend('Axial Direction', 'Radial Direction')
title({['Variance of linear fit for variable data points'],['Helix: k=' num2str(0)]})
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
title({['Variance of linear fit for variable data points'],['Helix: k=' num2str(0)]})
ylabel('Variance s^2 [A^2/m^2]')
xlabel('B_{app} [T] (Points Used)')



%% chi_a Simulation (High Fields)
load('N_simulation_High.mat')
chi_a_sim_High = inv(N_simulation_High); 

%% chi_a Simulation (Low Fields)
load('N_simulation_Low.mat')
chi_a_sim_Low = inv(N_simulation_Low);

%% chi_a Analytical (High Fields)

Beta = d_helix/L;
dir = ('Axial');
N_ax = DemagCylinder(Beta,dir)
N_rad = (1-N_ax)/2;

N_an = diag([N_rad N_rad N_ax])
chi_m = 24;
chi_a_an_High = inv(eye(3)/chi_m + N_an)

%% chi_a Experimental (Low Fields)

B_ax = 0.004; %Read from the Variance analysis plot
i_axial = length(find(H_app_ax_pos*mu_0 <= B_ax));

B_rad = 0.01; %Read from the Variance analysis plot
i_rad = length(find(H_app_rad_pos*mu_0 <= B_rad));

chi_a_Low_ax_ax = m_ax(i_axial,1);
chi_a_Low_ax_rad = m_ax(i_axial,2);
chi_a_Low_rad_ax = m_rad(i_rad,1);
chi_a_Low_rad_rad = m_rad(i_rad,2);

chi_a_Low = [...
    chi_a_Low_rad_rad   0                   chi_a_Low_ax_rad;...
    0                   chi_a_Low_rad_rad   chi_a_Low_ax_rad;...
    chi_a_Low_rad_ax    chi_a_Low_rad_ax    chi_a_Low_ax_ax];

%% Plot of Curves with Fits
subplot(2,1,1)
plot(H_app_ax_pos*mu_0,M_ax_ax_pos,'b-*',... % Magnetization ax,ax
    H_app_ax_pos*mu_0,M_ax_rad_pos,'r-*',... % Magnetization ax,rad
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_Low(3,3),'b--',... %Linearization (Low) ax,ax 
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_Low(2,3),'r--',...%Linearization (Low) ax,rad
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_sim_Low(3,3),'b-.',... %Simulation (Low) ax,ax 
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_sim_Low(2,3),...%Simulation (Low) ax,rad
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_sim_High(3,3),... %Simulation (High) ax,ax 
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_an_High(3,3),'k--',... %Analytical (High) ax,ax
    H_app_ax_pos*mu_0,ones(length(H_app_ax_pos),1)*M_s_ax_ax,'k-') %Saturation
    
axis tight
grid on
ylim([0 1.1*max([M_ax_ax_pos; M_ax_rad_pos])])
xlim([0 0.2])
legend('Axial Magnetization','Radial Magnetization','Axial Linearization','Radial Linearization', 'Axial Simulation (Low)','Radial Simulation (Low)','Axial Simulation (High)','Axial Analytical (High)','Axial Saturation')
title({'','Analysis of VSM Measurements (Wire)','','Simplified Magnetization Curves (Axial Direction)'})
xlabel('Applied Field B [T]')
ylabel('Magnetization M [A/m]')


subplot(2,1,2)
plot(H_app_rad_pos*mu_0,M_rad_rad_pos,'b-*',... % Magnetization rad,rad
    H_app_rad_pos*mu_0,M_rad_ax_pos,'r-*',... % Magnetization rad,ax
    H_app_rad_pos*mu_0,H_app_rad_pos*chi_a_Low(1,1),'b--',... %Linearization (Low) rad,rad 
    H_app_rad_pos*mu_0,H_app_rad_pos*chi_a_Low(3,1),'r--',...%Linearization (Low) rad,ax
    H_app_rad_pos*mu_0,H_app_rad_pos*chi_a_sim_Low(1,1),'b-.',... %Simulation (Low) rad,rad 
    H_app_rad_pos*mu_0,H_app_rad_pos*chi_a_sim_Low(3,1),...%Simulation (Low) rad,ax
    H_app_rad_pos*mu_0,H_app_rad_pos*chi_a_sim_High(1,1),... %Simulation (High) rad,rad 
    H_app_rad_pos*mu_0,H_app_rad_pos*chi_a_an_High(3,1),'k--',... %Analytical (High) rad,ax
    H_app_rad_pos*mu_0,ones(length(H_app_rad_pos),1)*M_s_rad_rad,'k-') %Saturation
    
axis tight
grid on
ylim([0 1.1*max([M_ax_ax_pos; M_ax_rad_pos])])
legend('Radial Magnetization','Axial Magnetization','Radial Linearization','Axial Linearization', 'Radial Simulation (Low)','Axial Simulation (Low)','Radial Simulation (High)','Raidal Analytical (High)','Radial Saturation')
title({'','Analysis of VSM Measurements (Wire)','','Simplified Magnetization Curves (Axial Direction)'})
xlabel('Applied Field B [T]')
ylabel('Magnetization M [A/m]')

%% Plot of Curves with Fits
subplot(2,1,1)
plot(H_app_ax_pos*mu_0,M_ax_ax_pos,'b-*',... % Magnetization ax,ax
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_Low(3,3),'b--',... %Linearization (Low) ax,ax 
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_sim_Low(3,3),'b-.',... %Simulation (Low) ax,ax
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_sim_High(3,3),'b:',... %Simulation (High) ax,ax 
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_an_High(3,3),'k--',... %Analytical (High) ax,ax
    H_app_ax_pos*mu_0,ones(length(H_app_ax_pos),1)*M_s_ax_ax,'k-') %Saturation
    
axis tight
grid on
ylim([0 1.1*max([M_ax_ax_pos; M_ax_rad_pos])])
xlim([0 0.2])
legend('Axial Magnetization','Axial Linearization', 'Axial Simulation (Low)','Axial Simulation (High)','Axial Analytical (High)','Axial Saturation')
title({'','Analysis of VSM Measurements (Wire)','','Simplified Magnetization Curves (Axial Direction)'})
xlabel('Applied Field B [T]')
ylabel('Magnetization M [A/m]')


subplot(2,1,2)
plot(H_app_ax_pos*mu_0,M_rad_rad_pos,'b-*',... % Magnetization rad,rad
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_Low(1,1),'b--',... %Linearization (Low) rad,rad 
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_sim_Low(1,1),'b-.',... %Simulation (Low) rad,rad 
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_sim_High(1,1),'b:',... %Simulation (High) rad,rad 
    H_app_ax_pos*mu_0,H_app_ax_pos*chi_a_an_High(1,1),'k--',... %Analytical (High) rad,rad
    H_app_ax_pos*mu_0,ones(length(H_app_ax_pos),1)*M_s_rad_rad,'k-') %Saturation
    
axis tight
grid on
ylim([0 1.1*max([M_ax_ax_pos; M_ax_rad_pos])])
legend('Radial Magnetization','Radial Linearization', 'Radial Simulation (Low)','Radial Simulation (High)','Radial Analytical (High)','Radial Saturation')
title({'','Analysis of VSM Measurements (Wire)','','Simplified Magnetization Curves (Axial Direction)'})
xlabel('Applied Field B [T]')
ylabel('Magnetization M [A/m]')


%% High Fields Saturation Ratio

% Axial Axial

H_app_ax_sat = M_s/chi_a_an_High(3,3);

M_ax_s_exp = interp1(H_app_ax_pos,M_ax_ax_pos,H_app_ax_sat,'spline');
r_ax = M_ax_s_exp/M_s

plot(H_app_ax_pos,M_ax_ax_pos,H_app_ax_pos,chi_a_an_High(3,3)*H_app_ax_pos,H_app_ax_pos,ones(length(H_app_ax_pos),1)*M_s)
hold on
scatter(H_app_ax_sat,M_ax_s_exp)
hold off
grid on
%ylim([min(M_ax_ax_pos) 1.1*M_ax_s_exp])
%xlim([min(H_app_ax_pos) 1.1*H_app_ax_sat])

%% Radial Radial

H_app_rad_sat = M_s/chi_a_an_High(1,1);

M_rad_s_exp = interp1(H_app_rad_pos,M_rad_rad_pos,H_app_rad_sat,'spline');
r_rad = M_rad_s_exp/M_s

plot(H_app_rad_pos,M_rad_rad_pos,H_app_rad_pos,chi_a_an_High(1,1)*H_app_rad_pos,H_app_rad_pos,ones(length(H_app_rad_pos),1)*M_s)
hold on
scatter(H_app_rad_sat,M_rad_s_exp)
hold off
grid on
ylim([min(M_rad_rad_pos) 1.1*M_rad_s_exp])
%xlim([min(H_app_rad_pos) 1.1*H_app_rad_sat])