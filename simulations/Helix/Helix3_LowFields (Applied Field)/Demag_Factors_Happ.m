close all
clear all
clc

Data_M_x = dlmread('M_x.txt');
Data_M_y = dlmread('M_y.txt');
Data_M_z = dlmread('M_z.txt');
Data_H_x = dlmread('H_x.txt');
Data_H_y = dlmread('H_y.txt');
Data_H_z = dlmread('H_z.txt');
Data_B_app = dlmread('B_app.txt');
Data_n_helix = dlmread('n_helix.txt');
mu_0 = 1.25663706e-6;
mu_r = 1;
mu = mu_r*mu_0;

n_helix =Data_n_helix(1,5);

[row, col] = size(Data_B_app);
n_max = row/3;

%%

for i = 1 : n_max
B_app = Data_B_app(3*(i-1)+1,1);
H_app_val = B_app/mu; %T
B_app_vec(i) = B_app;

%B_x
M_x_x = Data_M_x(3*(i-1)+1,5); %A/m
M_x_y = Data_M_y(3*(i-1)+1,5); %A/m
M_x_z = Data_M_z(3*(i-1)+1,5); %A/m

M_x = [M_x_x M_x_y M_x_z]';

H_x_x = Data_H_x(3*(i-1)+1,5); %A/m
H_x_y = Data_H_y(3*(i-1)+1,5);  %A/m
H_x_z = Data_H_z(3*(i-1)+1,5);  %A/m

H_x = [H_x_x H_x_y H_x_z]';
H_app_x = [H_app_val 0 0]';

%B_y
M_y_x = Data_M_x(3*(i-1)+2,5); %A/m
M_y_y = Data_M_y(3*(i-1)+2,5); %A/m
M_y_z = Data_M_z(3*(i-1)+2,5); %A/m

M_y = [M_y_x M_y_y M_y_z]';

H_y_x = Data_H_x(3*(i-1)+2,5); %A/m
H_y_y = Data_H_y(3*(i-1)+2,5);  %A/m
H_y_z = Data_H_z(3*(i-1)+2,5);  %A/m

H_y = [H_y_x H_y_y H_y_z]';
H_app_y = [0 H_app_val 0]';


%B_z
M_z_x = Data_M_x(3*(i-1)+3,5); %A/m
M_z_y = Data_M_y(3*(i-1)+3,5); %A/m
M_z_z = Data_M_z(3*(i-1)+3,5); %A/m

M_z = [M_z_x M_z_y M_z_z]';


H_z_x = Data_H_x(3*(i-1)+3,5); %A/m
H_z_y = Data_H_y(3*(i-1)+3,5); %A/m
H_z_z = Data_H_z(3*(i-1)+3,5); %A/m

H_z = [H_z_x H_z_y H_z_z]';

H_app_z = [0 0 H_app_val]';


%Demag Factors
M = [M_x M_y M_z];
H = [H_x H_y H_z];
H_app = [H_app_x H_app_y H_app_z];


N = (H_app - H)*inv(M);
TR_vec(i) = trace(N);

[V,D] = eig(N);


N_11_vec(i) = N(1,1);
N_22_vec(i) = N(2,2);
N_33_vec(i) = N(3,3);

Results.(['H' num2str(n_helix)]).(['B' num2str(B_app*1000) 'mT']).N = N;


end

figure(3)
subplot(4,1,1)
plot(B_app_vec, N_11_vec)
title('N_{11}')
grid on
xlabel('B_{app} [T]')

subplot(4,1,2)
plot(B_app_vec, N_22_vec)
title('N_{22}')
grid on
xlabel('Helix')
xlabel('B_{app} [T]')

subplot(4,1,3)
plot(B_app_vec, N_33_vec)
title('N_{33}')
grid on
xlabel('B_{app} [T]')

subplot(4,1,4)
plot(B_app_vec, TR_vec)
title('tr(N)')
grid on
xlabel('B_{app} [T]')

%% Saving

save('Results.mat','Results')


