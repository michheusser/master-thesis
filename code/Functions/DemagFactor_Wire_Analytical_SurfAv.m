function [N_matrix_global_av ,N_matrix_local_av, rot_Matrix_local] = DemagFactor_Wire_Analytical_SurfAv(s, N_R, N_phi, shape, shape_param)
% Calculates the surface average N-Tensor at a certain point

N_matrix_global_av = zeros(3,3);
N_matrix_local_av = zeros(3,3);

phi_vec = [0:N_phi-1]*2*pi/N_phi;

if(strcmp(shape,'thin helix'))
    n_helix = shape_param;
    [~, ~, ~, d_helix] = Helixinfo(n_helix);
    thickness = 3e-7
    R_min = d_helix/2 - thickness;
    R_max = d_helix/2;
elseif(strcmp(shape,'helix'))
    n_helix = shape_param;
    [~, ~, ~, d_helix] = Helixinfo(n_helix);
    R_min = 0;
    R_max = d_helix/2;
elseif(strcmp(shape,'rod'))
    R_wire = shape_param(1);
    R_min = 0;
    R_max = R_wire;
end

R_vec = R_min + [1: N_R]*(R_max - R_min)/(N_R+1);

disp(['Surface Average at s = ' num2str(s)])
for i_phi = 1 : N_phi
    for i_R = 1 : N_R
        phi_i = phi_vec(i_phi);
        R_i = R_vec(i_R);
        s_i = s;
        
        coord = [phi_i, R_i, s_i];
        [N_matrix_global,N_matrix_local, rot_Matrix_local] = DemagFactor_Wire_Analytical(coord,shape,shape_param);
        
        N_matrix_global_av = N_matrix_global_av + 1/(N_R*N_phi)*N_matrix_global;
        N_matrix_local_av = N_matrix_local_av + 1/(N_R*N_phi)*N_matrix_local;
        disp(['Surface Average Progress: ' num2str(round(((i_phi-1)*N_R+i_R)*100/(N_phi*N_R)) ) '%'])
    end
    
end
