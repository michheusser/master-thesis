function rot_Matrix_local = Rot_Matrix(s_i,shape,shape_param)

% Calculates the local rotation matrix at a specific point inside a helix
% of helix index n_helix

% shape: 'helix', 'rod', 'other'
% Helix: shape_param = n_helix
% Thin Helix: shape_param = n_helix
% Rod: shape_param = [R_wire, L_rod]
% Other: shape_param = [f_fun(s_param), R_wire, s_min, s_max, R_rot(s_param)]

syms s_param s s_str phi phi_str R R_str


if(strcmp(shape,'helix'))
    
    n_helix = shape_param;
    [h_helix, D_helix,~, ~] = Helixinfo(n_helix);
    f_fun = [D_helix/2 * cos(s_param); D_helix/2 * sin(s_param); h_helix*s_param/(2*pi)];
    dfds = diff(f_fun,s_param);
    z_w = dfds/norm(dfds);
    x_w = [cos(s_param) -sin(s_param) 0; sin(s_param) cos(s_param) 0; 0 0 1]*[1; 0; 0];
    y_w = cross(z_w,x_w);
    rot_Matrix(s_param) = [x_w y_w z_w];
    
    
elseif(strcmp(shape,'thin helix'))
    
    n_helix = shape_param;
    [h_helix, D_helix, k_helix, d_helix] = Helixinfo(n_helix);
    R_wire = d_helix/2;
    s_min = 0;
    s_max = k_helix*2*pi;
    f_fun = [D_helix/2 * cos(s_param); D_helix/2 * sin(s_param); h_helix*s_param/(2*pi)];
    dfds = diff(f_fun,s_param);
    z_w = dfds/norm(dfds);
    x_w = [cos(s_param) -sin(s_param) 0; sin(s_param) cos(s_param) 0; 0 0 1]*[1; 0; 0];
    y_w = cross(z_w,x_w);
    rot_Matrix(s_param) = [x_w y_w z_w];
    
    thickness = 3e-7;
    R_min = R_wire - thickness;
    R_max = R_wire;
    
elseif(strcmp(shape,'rod'))
    R_wire = shape_param(1);
    L_rod = shape_param(2);
    s_min = 0;
    s_max = L_rod;
    f_fun = [0; 0; s_param];
    dfds = diff(f_fun,s_param);
    z_w = dfds/norm(dfds);
    x_w = [1 0 0]';
    y_w = [0 1 0]';
    rot_Matrix(s_param) = [x_w y_w z_w];
    
end



rot_Matrix_local = double(rot_Matrix(s_i));
%N_matrix_local = rot_Matrix_local'*N_matrix_global*rot_Matrix_local;

