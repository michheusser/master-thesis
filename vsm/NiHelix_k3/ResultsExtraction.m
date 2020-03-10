close all
clear all
clc

%%
k_helix = 3;

mu_0 = 1.25663706e-6;
Oe_to_Am = 1e3/(4*pi);
emu_to_Am2 = 1e-3;

d_helix = 250e-6;
D_helix = 1.51e-3 + d_helix;
h_helix = 0.4e-3;
L = k_helix*sqrt(h_helix^2+(D_helix*pi)^2);
A = (d_helix/2)^2*pi;
V = L*A;

Data_ax = dlmread('Clean/Helix_Ni_k3-VIR_axial.txt', ',', 1, 0);
H_app_ax = Data_ax(:,8)*Oe_to_Am;
M_ax_ax = (Data_ax(:,10)/1000)*emu_to_Am2/V;
M_ax_rad = (Data_ax(:,11)/1000)*emu_to_Am2/V;

Data_rad = dlmread('Clean/Helix_Ni_k3-VIR_radial.txt', ',', 1, 0);
H_app_rad = Data_rad(:,8)*Oe_to_Am;
M_rad_rad = (Data_rad(:,10)/1000)*emu_to_Am2/V;
M_rad_ax = -(Data_rad(:,11)/1000)*emu_to_Am2/V;

plot(H_app_ax*mu_0,M_ax_ax,H_app_ax*mu_0,M_ax_rad,H_app_rad*mu_0,M_rad_ax,H_app_rad*mu_0,M_rad_rad)
grid on
axis tight
xlabel('Applied Field B [T]')
ylabel('Magnetization M [A/m]')
legend('M_{ax,ax}','M_{ax,rad}','M_{rad,ax}','M_{rad,rad}')

%% Removal of initial magnetization
n_start = 100
n_end = 491

H_app_ax_loop = H_app_ax(n_start:n_end); 
M_ax_ax_loop = M_ax_ax(n_start:n_end);
M_ax_rad_loop = M_ax_rad(n_start:n_end);

H_app_rad_loop = H_app_rad(n_start:n_end);
M_rad_rad_loop = M_rad_rad(n_start:n_end);
M_rad_ax_loop = M_rad_ax(n_start:n_end);

plot(H_app_ax_loop*mu_0,M_ax_ax_loop,H_app_ax_loop*mu_0,M_ax_rad_loop,H_app_rad_loop*mu_0,M_rad_ax_loop,H_app_rad_loop*mu_0,M_rad_rad_loop)
grid on
axis tight
xlabel('Applied Field B [T]')
ylabel('Magnetization M [A/m]')
legend('M_{ax,ax}','M_{ax,rad}','M_{rad,ax}','M_{rad,rad}')

%% Low field truncation
B_app_max = 0.005;
H_app_max = B_app_max/mu_0;

H_app_max_ax_ax = H_app_max;
H_app_max_ax_rad = H_app_max;
H_app_max_rad_ax = H_app_max;
H_app_max_rad_rad = H_app_max;

I_low_ax_ax = find( (abs(H_app_ax_loop) <= H_app_max_ax_ax));
I_low_ax_rad = find( (abs(H_app_ax_loop) <= H_app_max_ax_rad));
I_low_rad_ax = find( (abs(H_app_rad_loop) <= H_app_max_rad_ax));
I_low_rad_rad = find( (abs(H_app_rad_loop) <= H_app_max_rad_rad));

H_app_ax_ax_low = H_app_ax_loop(I_low_ax_ax); 
M_ax_ax_low = M_ax_ax_loop(I_low_ax_ax);

H_app_ax_rad_low = H_app_ax_loop(I_low_ax_rad); 
M_ax_rad_low = M_ax_rad_loop(I_low_ax_rad);

H_app_rad_ax_low = H_app_rad_loop(I_low_rad_ax);
M_rad_ax_low = M_rad_ax_loop(I_low_rad_ax);

H_app_rad_rad_low = H_app_rad_loop(I_low_rad_rad);
M_rad_rad_low = M_rad_rad_loop(I_low_rad_rad);

M_ax_ax_fit = fit(H_app_ax_ax_low,M_ax_ax_low,'poly1');
M_ax_rad_fit = fit(H_app_ax_rad_low,M_ax_rad_low,'poly1');
M_rad_ax_fit = fit(H_app_rad_ax_low,M_rad_ax_low,'poly1');
M_rad_rad_fit = fit(H_app_rad_rad_low,M_rad_rad_low,'poly1');


plot(...
    H_app_ax_ax_low*mu_0,M_ax_ax_low,'b--*',...
    H_app_ax_ax_low*mu_0,M_ax_ax_fit(H_app_ax_ax_low),'b-',...
    H_app_ax_rad_low*mu_0,M_ax_rad_low,'g--*',...
    H_app_ax_rad_low*mu_0,M_ax_rad_fit(H_app_ax_rad_low),'g-',...
    H_app_rad_ax_low*mu_0,M_rad_ax_low,'r--*',...
    H_app_rad_ax_low*mu_0,M_rad_ax_fit(H_app_rad_ax_low),'r-',...
    H_app_rad_rad_low*mu_0,M_rad_rad_low,'k--*',...
    H_app_rad_rad_low*mu_0,M_rad_rad_fit(H_app_rad_rad_low),'k-')

grid on
axis tight
xlabel('Applied Field B [T]')
ylabel('Magnetization M [A/m]')
legend('M_{ax,ax}','M_{ax,ax: Fitted}','M_{ax,rad}','M_{ax,rad: Fitted}','M_{rad,ax}','M_{rad,ax: Fitted}','M_{rad,rad}','M_{rad,rad: Fitted}')

% chi_a
M_ax_ax_low_coeff = coeffvalues(M_ax_ax_fit);
chi_a_ax_ax = M_ax_ax_low_coeff(1);

M_ax_rad_low_coeff = coeffvalues(M_ax_rad_fit);
chi_a_ax_rad = M_ax_rad_low_coeff(1);

M_rad_ax_low_coeff = coeffvalues(M_rad_ax_fit);
chi_a_rad_ax = M_rad_ax_low_coeff(1);

M_rad_rad_low_coeff = coeffvalues(M_rad_rad_fit);
chi_a_rad_rad = M_rad_rad_low_coeff(1);

chi_a = [...
    chi_a_ax_ax chi_a_rad_ax chi_a_rad_ax;...
    chi_a_ax_rad chi_a_rad_rad 0;...
    chi_a_ax_rad 0 chi_a_rad_rad];

N_matrix = inv(chi_a);
N_matrix_global = circshift(N_matrix,[-1 -1]);

N_simulationB
Trace = trace(N_simulationB)
r = N_simulationB(3,3)/N_simulationB(1,1)
N_matrix_global
Trace = trace(N_matrix_global)
r = N_matrix_global(3,3)/N_matrix_global(1,1)


%% Simulation

mphstart()
%%
model = mphload('Helix_VSM_k3_Nickel.mph')
%%
mu_0 = 1.25663706e-6;
mu_r = 1;
mu = mu_r*mu_0;


    [H_x, H_y, H_z, M_x, M_y, M_z, B_app] = mphmean(model,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'B_app'},'volume','selection',2,'dataset','dset1');
    
    H = [H_x' H_y' H_z'];
    M = [M_x' M_y' M_z'];
    H_app = B_app(1)/mu*eye(3);
    
    N_simulationB = (H_app - H)/M;
    DemagFactors.N_simulationB = N_simulationB;
    DemagFactors.N_matrix = N_matrix_global;


save('DemagFactors.mat','DemagFactors')

%%
N_matrix_global
N_simulationB
