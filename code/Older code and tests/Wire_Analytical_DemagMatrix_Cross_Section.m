close all
clear all
clc


% %HELIX
% n_helix = 3;
% shape = 'helix';
% shape_param = n_helix;
% n = 10;
% %r_i_vec = [0:n-1]/n*1.5e-6;
% r_i_vec = 0;
% phi_i = 0;
% %s_i_vec = [1:n-1]'/n*3*2*pi;
% s_i = 3*2*pi/2;

%THIN HELIX
n_helix = 3;
shape = 'thin helix';
shape_param = n_helix;
thickness = 0.3e-6; %300 nm
R_wire = 1.5e-6;
r_i = R_wire - 0.5*thickness;
phi_i = 0;
n = 20;
s_i_vec = [1:n-1]'/n*3*2*pi;

% %ROD
% R_wire = 0.1;
% L_rod = 1;
% r_i_vec = 0;
% r_i = R_wire;
% phi_i = 0;
% n = 20;
% %s_i_vec = [0:n]'/(n)*L_rod;
% s_i_vec = L_rod/2;
% shape = 'rod';
% shape_param = [R_wire, L_rod];


N = length(s_i_vec);

N_11_local_vec = zeros(N,1);
N_22_local_vec = zeros(N,1);
N_33_local_vec = zeros(N,1);

N_11_global_vec = zeros(N,1);
N_22_global_vec = zeros(N,1);
N_33_global_vec = zeros(N,1);

N_local_cell = cell(N,1);
N_global_cell = cell(N,1);

N_matrix_global_average = 0;
%for i = 1 : length(r_i_vec)
for i = 1 : length(s_i_vec)
    i
    %r_i = r_i_vec(i);
    s_i = s_i_vec(i);
    coord = [phi_i r_i s_i];
    
    [N_matrix_global, N_matrix_local, rot_Matrix_local] = DemagFactor_Wire_Analytical(coord,shape,shape_param);
    
    N_local_cell{i,1} = N_matrix_local;
    N_11_local_vec(i) = N_matrix_local(1,1);
    N_22_local_vec(i) = N_matrix_local(2,2);
    N_33_local_vec(i) = N_matrix_local(3,3);
    
    N_matrix_local
    TR = trace(N_matrix_local)
    
    N_global_cell{i,1} = N_matrix_local;
    N_11_global_vec(i) = N_matrix_global(1,1);
    N_22_global_vec(i) = N_matrix_global(2,2);
    N_33_global_vec(i) = N_matrix_global(3,3);
   
    N_matrix_global_average = N_matrix_global_average + N_matrix_global/length(s_i_vec);
    
end



%%
vec_to_plot = s_i_vec;


figure(1)
plot(vec_to_plot,N_11_local_vec,'*-',vec_to_plot,N_22_local_vec,'o-',vec_to_plot, N_33_local_vec,'x-')
grid on
axis tight
xlabel('s')
ylabel('N')
%ylim([0 1])
title(['Demagnetization Factors along the cross section of a ' shape ' at the bottom (Local Coordinates)'])
legend('N_{helix radial}', 'N_{y}','N_{wire axis}','Location','Best')

figure(2)
plot(vec_to_plot,N_11_global_vec,'*-',vec_to_plot,N_22_global_vec,'o-',vec_to_plot, N_33_global_vec,'x-')
grid on
axis tight
xlabel('s')
ylabel('s')
%ylim([0 1])
title(['Demagnetization Factors along the cross section of a ' shape ' at the bottom (Global Coordinates)'])
legend('N_{x}', 'N_{y}','N_{z}','Location','Best')


