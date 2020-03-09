close all
clear all
clc

n_helix_curr = 9;
n_helix_vec_long = [1:0.1:3 4:10];

n_helix_vec = n_helix_vec_long(n_helix_curr);

N_s = 20;
N_R = 3;
N_phi = 5;

n_cores = 12;

matlabpool('local',n_cores)

for j = 1 : length(n_helix_vec)    

n_helix = n_helix_vec(j);


%THIN HELIX
shape = 'thin helix';
shape_param = n_helix;
n = N_s;
s_i_vec = [1:n]'/(n+1)*3*2*pi;

N_matrix_global_surf_av_cell = cell(length(s_i_vec),1);
parfor i = 1 : length(s_i_vec)
    s_i = s_i_vec(i);
    [N_matrix_global_surf_av ,N_matrix_local_surf_av, rot_Matrix_local] = DemagFactor_Wire_Analytical_SurfAv(s_i, N_R, N_phi, shape, shape_param);
    N_matrix_global_surf_av_cell{i,1} = N_matrix_global_surf_av;
end

N_matrix_global_vol_av = 0;
for i = 1 : length(s_i_vec)
   
    N_matrix_global_vol_av = N_matrix_global_vol_av + N_matrix_global_surf_av_cell{i,1}/length(s_i_vec);
    
end

DemagMatrix_Vol_Average_Analytical.(['S' num2str(j)]).n_helix = n_helix;
DemagMatrix_Vol_Average_Analytical.(['S' num2str(j)]).N_matrix_global_vol_cell = N_matrix_global_surf_av_cell;
DemagMatrix_Vol_Average_Analytical.(['S' num2str(j)]).N_matrix_global_vol_av = N_matrix_global_vol_av;




end


save(['DemagMatrix_ThinHelix_Vol_Average_Analytical_' num2str(n_helix_curr) '.mat'],'DemagMatrix_Vol_Average_Analytical');
matlabpool close

