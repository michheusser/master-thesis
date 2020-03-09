close all
clear all
clc

%
shape = 'helix';
n_helix = 9;
    shape = 'helix';
    shape_param = n_helix;
r_i = 0;
phi_i = 0;
n = 20;
s_i_vec = [1:n-1]'/n*3*2*pi;

N_11_vec = zeros(length(s_i_vec),1);
N_22_vec = zeros(length(s_i_vec),1);
N_33_vec = zeros(length(s_i_vec),1);

N_cell = cell(length(s_i_vec),1);
N_cyl_cell = cell(length(s_i_vec),1);

for i = 1 : length(s_i_vec)
    i
    s_i = s_i_vec(i);
    coord = [phi_i r_i s_i];

    [N_matrix_local, rot_Matrix_local] = DemagFactor_Wire_Analytical(coord,shape,shape_param);
    
    N_cell{i,1} = N_matrix_local;
    N_11_vec(i) = N_matrix_local(1,1);
    N_22_vec(i) = N_matrix_local(2,2);
    N_33_vec(i) = N_matrix_local(3,3);
    
end

%%

figure(1)
plot(s_i_vec/(2*pi),N_11_vec,'*-',s_i_vec/(2*pi),N_22_vec,'o-',s_i_vec/(2*pi), N_33_vec,'x-')
grid on
xlabel('Helix')
ylabel('N')
%ylim([0 1])
title('Demagnetization Factors along the center of a helix')
legend('N_{radial}', 'N_{2}','N_{tangential}','Location','Best')




