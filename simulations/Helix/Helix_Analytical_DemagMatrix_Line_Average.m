close all
clear all
clc

n_helix_vec = [1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 4 5 6 7 8 9];

N_s = 20;

n_cores = 12;

matlabpool('local',n_cores)

for j = 1 : length(n_helix_vec)    

n_helix = n_helix_vec(j);


%HELIX
shape = 'helix';
shape_param = n_helix;
n = N_s;
s_i_vec = [1:n]'/(n+1)*3*2*pi;

N_matrix_global_line_cell = cell(length(s_i_vec),1);
parfor i = 1 : length(s_i_vec)
    s_i = s_i_vec(i);
    phi_i = 0;
    r_i = 0;
    coord = [phi_i r_i s_i];
    [N_matrix_global ,N_matrix_local, rot_Matrix_local] = DemagFactor_Wire_Analytical(coord, shape, shape_param);
    N_matrix_global_line_cell{i,1} = N_matrix_global;
end

N_matrix_global_line_av = 0;
for i = 1 : length(s_i_vec)
   
    N_matrix_global_line_av = N_matrix_global_line_av + N_matrix_global_line_cell{i,1}/length(s_i_vec);
    
end

DemagMatrix_Line_Average_Analytical.(['S' num2str(j)]).n_helix = n_helix;
DemagMatrix_Line_Average_Analytical.(['S' num2str(j)]).N_matrix_global_line_cell = N_matrix_global_line_cell;
DemagMatrix_Line_Average_Analytical.(['S' num2str(j)]).N_matrix_global_line_av = N_matrix_global_line_av;




end


%save('DemagMatrix_Helix_Line_Average_Analytical.mat','DemagMatrix_Line_Average_Analytical');
matlabpool close

