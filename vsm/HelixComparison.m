close all
clear all
clc

%% Data Extraction

k_vec = [2 3 4 7];

% k = 2
N_simulation_High_data2 = load('NiHelix_k2/N_simulation_High.mat')
N_simulation_High2 = N_simulation_High_data2.N_simulation_High;

N_simulation_Low_data2 = load('NiHelix_k2/N_simulation_Low.mat')
N_simulation_Low2 = N_simulation_Low_data2.N_simulation_Low;

N_VSM_Low_data2 = load('NiHelix_k2/N_VSM_Low.mat')
N_VSM_Low2 = N_VSM_Low_data2.N_VSM_Low;

% k = 3
N_simulation_High_data3 = load('NiHelix_k3/N_simulation_High.mat')
N_simulation_High3 = N_simulation_High_data3.N_simulation_High;

N_simulation_Low_data3 = load('NiHelix_k3/N_simulation_Low.mat')
N_simulation_Low3 = N_simulation_Low_data3.N_simulation_Low;

N_VSM_Low_data3 = load('NiHelix_k3/N_VSM_Low.mat')
N_VSM_Low3 = N_VSM_Low_data3.N_VSM_Low;

% k = 4
N_simulation_High_data4 = load('NiHelix_k4/N_simulation_High.mat')
N_simulation_High4 = N_simulation_High_data4.N_simulation_High;

N_simulation_Low_data4 = load('NiHelix_k4/N_simulation_Low.mat')
N_simulation_Low4 = N_simulation_Low_data4.N_simulation_Low;

N_VSM_Low_data4 = load('NiHelix_k4/N_VSM_Low.mat')
N_VSM_Low4 = N_VSM_Low_data4.N_VSM_Low;

% k = 7
N_simulation_High_data7 = load('NiHelix_k7/N_simulation_High.mat')
N_simulation_High7 = N_simulation_High_data7.N_simulation_High;

N_simulation_Low_data7 = load('NiHelix_k7/N_simulation_Low.mat')
N_simulation_Low7 = N_simulation_Low_data7.N_simulation_Low;

N_VSM_Low_data7 = load('NiHelix_k7/N_VSM_Low.mat')
N_VSM_Low7 = N_VSM_Low_data7.N_VSM_Low;

%%
N_11_sim_High = zeros(4,1);
N_22_sim_High = zeros(4,1);
N_33_sim_High = zeros(4,1);
N_11_sim_Low = zeros(4,1);
N_22_sim_Low = zeros(4,1);
N_33_sim_Low = zeros(4,1);

N_11_VSM_Low = zeros(4,1);
N_22_VSM_Low = zeros(4,1);
N_33_VSM_Low = zeros(4,1);

N_11_sim_High(1) = N_simulation_High2(1,1);
N_11_sim_High(2) = N_simulation_High3(1,1);
N_11_sim_High(3) = N_simulation_High4(1,1);
N_11_sim_High(4) = N_simulation_High7(1,1);

N_22_sim_High(1) = N_simulation_High2(2,2);
N_22_sim_High(2) = N_simulation_High3(2,2);
N_22_sim_High(3) = N_simulation_High4(2,2);
N_22_sim_High(4) = N_simulation_High7(2,2);

N_33_sim_High(1) = N_simulation_High2(3,3);
N_33_sim_High(2) = N_simulation_High3(3,3);
N_33_sim_High(3) = N_simulation_High4(3,3);
N_33_sim_High(4) = N_simulation_High7(3,3);

N_11_sim_Low(1) = N_simulation_Low2(1,1);
N_11_sim_Low(2) = N_simulation_Low3(1,1);
N_11_sim_Low(3) = N_simulation_Low4(1,1);
N_11_sim_Low(4) = N_simulation_Low7(1,1);

N_22_sim_Low(1) = N_simulation_Low2(2,2);
N_22_sim_Low(2) = N_simulation_Low3(2,2);
N_22_sim_Low(3) = N_simulation_Low4(2,2);
N_22_sim_Low(4) = N_simulation_Low7(2,2);

N_33_sim_Low(1) = N_simulation_Low2(3,3);
N_33_sim_Low(2) = N_simulation_Low3(3,3);
N_33_sim_Low(3) = N_simulation_Low4(3,3);
N_33_sim_Low(4) = N_simulation_Low7(3,3);

N_11_VSM_Low(1) = N_VSM_Low2(1,1);
N_11_VSM_Low(2) = N_VSM_Low3(1,1);
N_11_VSM_Low(3) = N_VSM_Low4(1,1);
N_11_VSM_Low(4) = N_VSM_Low7(1,1);

N_22_VSM_Low(1) = N_VSM_Low2(2,2);
N_22_VSM_Low(2) = N_VSM_Low3(2,2);
N_22_VSM_Low(3) = N_VSM_Low4(2,2);
N_22_VSM_Low(4) = N_VSM_Low7(2,2);

N_33_VSM_Low(1) = N_VSM_Low2(3,3);
N_33_VSM_Low(2) = N_VSM_Low3(3,3);
N_33_VSM_Low(3) = N_VSM_Low4(3,3);
N_33_VSM_Low(4) = N_VSM_Low7(3,3);

subplot(3,1,1)
plot(k_vec,N_11_sim_High,'*--',k_vec,N_11_sim_Low,'*--',k_vec,N_11_VSM_Low,'*--')
title('Demagnetization Factors for Measured Helices')
xlabel('k (Number of Coils)')
ylabel('N_{11}')
legend('High Fields: Simulation','Low Fields: Simulation','Low Fields: VSM')
subplot(3,1,2)
plot(k_vec,N_22_sim_High,'*--',k_vec,N_22_sim_Low,'*--',k_vec,N_22_VSM_Low,'*--')
xlabel('k (Number of Coils)')
ylabel('N_{22}')
legend('High Fields: Simulation','Low Fields: Simulation','Low Fields: VSM')
subplot(3,1,3)
plot(k_vec,N_33_sim_High,'*--',k_vec,N_33_sim_Low,'*--',k_vec,N_33_VSM_Low,'*--')
ylabel('N_{33}')
xlabel('k (Number of Coils)')
legend('High Fields: Simulation','Low Fields: Simulation','Low Fields: VSM')