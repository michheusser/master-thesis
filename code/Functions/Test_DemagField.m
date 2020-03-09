close all
clear all
clc

coord = [0 0 0.5];
shape = 'rod';
shape_param = [0.1 1];

 factor = 0.0000;
% M_x = @(R_param,phi_param,s_param)(1+factor*R_param^2)*cos(1+factor*phi_param^2);
% M_y = @(R_param,phi_param,s_param)(1+factor*R_param^2)*sin(1+factor*phi_param^2);
% M_z = @(R_param,phi_param,s_param) 1+factor*s_param^2;

 M_x = @(R_param,phi_param,s_param) 1; 
 M_y = @(R_param,phi_param,s_param) 0; 
 M_z = @(R_param,phi_param,s_param) 0; 
 
H_demag_x = DemagField_Wire_Analytical(coord, shape, shape_param, M_x, M_y, M_z)

 M_x = @(R_param,phi_param,s_param) 0;
 M_y = @(R_param,phi_param,s_param) 1; 
 M_z = @(R_param,phi_param,s_param) 0;

H_demag_y = DemagField_Wire_Analytical(coord, shape, shape_param, M_x, M_y, M_z)


 M_x = @(R_param,phi_param,s_param) 0;
 M_y = @(R_param,phi_param,s_param) 0;
 M_z = @(R_param,phi_param,s_param) 1;

H_demag_z = DemagField_Wire_Analytical(coord, shape, shape_param, M_x, M_y, M_z)


N_approx = -[H_demag_x H_demag_y H_demag_z] 