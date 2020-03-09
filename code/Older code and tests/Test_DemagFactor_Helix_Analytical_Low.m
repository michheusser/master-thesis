close all
clear all
clc

mu_0 = 1.25663706e-6;
mu_r = 1;
mu = mu_r*mu_0;

N = 10;
n_helix = 3;
chi_m = 1000;
shape_function = 'fourier';
integration_method = 'auto';
[h_helix, D_helix, k_helix, d_helix] = Helixinfo(n_helix);

s_min = 0;
s_max = 2*pi*k_helix;

N_points = 50;
s_vec = s_min + [1:N_points]/(N_points+1)*(s_max-s_min);

[N_global, N_local, N_average, N_average_2] = DemagFactor_Helix_Analytical_Low_Numerical2(n_helix,chi_m,N,shape_function, integration_method, s_vec)

%% Simulation
mphstart()
model = mphload('Helix_LowFields (Applied Field)2.mph')

i_helix = 11;

n_helix_vec = [1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 4 5 6 7 8 9 10];
n_helix = n_helix_vec(i_helix);

syms s phi R
r_c(s) = [(D_helix/2)*cos(s); (D_helix/2)*sin(s); h_helix/(2*pi)*s]+[0 ; 0 ; -k_helix*h_helix/2];
dfds(s) = diff(r_c(s),s);
z_w(s) = dfds/norm(dfds(s));
x_w(s) = [cos(s) -sin(s) 0; sin(s) cos(s) 0; 0 0 1]*[1; 0; 0];
y_w(s) = cross(z_w(s),x_w(s));
rot_Matrix(s) = [x_w(s) y_w(s) z_w(s)];

r(R,phi,s) = rot_Matrix(s)*[R*cos(phi); R*sin(phi); 0]+r_c(s);
%%

R_i = 0;
phi_i = 0;
COORD = zeros(3,length(s_i_vec));

for i = 1 : length(s_i_vec)
    s_i = s_vec(i);
    COORD(:,i) = double(r(R_i,phi_i,s_i));
end


[H_x, H_y, H_z, M_x, M_y, M_z, B_app, n_helix_vec] = mphinterp(model,{'mfnc.Hx', 'mfnc.Hy', 'mfnc.Hz', 'mfnc.Mx', 'mfnc.My', 'mfnc.Mz', 'B_app', 'n_helix'},'coord',COORD,'selection',2,'dataset','dset2','outersolnum',i_helix);

N_simulation_global_cell = cell(length(s_i_vec),1);
N_simulation_local_cell = cell(length(s_i_vec),1);

for i = 1 : length(s_i_vec)
    i
    H_x_i = H_x(:,i);
    H_y_i = H_y(:,i);
    H_z_i = H_z(:,i);
    M_x_i = M_x(:,i);
    M_y_i = M_y(:,i);
    M_z_i = M_z(:,i);
    
    H = [H_x_i'; H_y_i'; H_z_i'];
    M = [M_x_i'; M_y_i'; M_z_i'];
    H_app = B_app(1,1)/mu*eye(3);
    
    N_simulationB = (H_app - H)*inv(M);
    N_simulation_global_cell{i,1} = N_simulationB;
    N_simulation_local_cell{i,1} = double(rot_Matrix(s_vec(i)).'*N_simulation_global_cell{i,1}*rot_Matrix(s_vec(i)));
    
end


%% Plot
close all
s_plot = s_vec;
N_an_11_vec = zeros(length(s_plot),1);
N_an_22_vec = zeros(length(s_plot),1);
N_an_33_vec = zeros(length(s_plot),1);

N_sim_11_vec = zeros(length(s_plot),1);
N_sim_22_vec = zeros(length(s_plot),1);
N_sim_33_vec = zeros(length(s_plot),1);

for i = 1 : length(s_plot)
    s_i = s_plot(i);
    N_an_11_vec(i) = double(N_local{i,1}(1,1));
    N_an_22_vec(i) = double(N_local{i,1}(2,2));
    N_an_33_vec(i) = double(N_local{i,1}(3,3));
    
    %N_sim_11_vec(i) = N_simulation_local_cell{i,1}(1,1);
    %N_sim_22_vec(i) = N_simulation_local_cell{i,1}(2,2);
    %N_sim_33_vec(i) = N_simulation_local_cell{i,1}(3,3);
    
end

plot(s_plot,N_an_11_vec,'Color','r','LineStyle','-')
hold on
plot(s_plot,N_sim_11_vec,'Color','r','LineStyle','-.')
plot(s_plot,N_an_22_vec,'Color','b','LineStyle','-')
plot(s_plot,N_sim_22_vec,'Color','b','LineStyle','-.')
plot(s_plot,N_an_33_vec,'Color','g','LineStyle','-')
plot(s_plot,N_sim_33_vec,'Color','g','LineStyle','-.')
hold off
legend('N_{analytical,11}(s)','N_{simulation,11}(s)','N_{analytical,22}(s)','N_{simulation,22}(s)','N_{analytical33}(s)','N_{simulation,33}(s)')
grid on
axis tight
xlabel('s')

title({'',['Demagnetization Factor for Low Fields (N = ' num2str(N) ')'],...
    ['n_{helix} = ' num2str(n_helix)],['Method: ' shape_function ' (N = ' num2str(N) ')'],['\chi_m = ' num2str(chi_m)]})
