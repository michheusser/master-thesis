%function N_matrix = DemagFactor_Helix_Analytical(n_helix,s_helix)
% Calculates the Demagnetization Matrix at a specific point inside a helix
% of helix index n_helix

%Source: "A Sum Rule Concerning the Inhomogeneous Demagnetizing Field in
%Nonellipsoidal Samples, Ernst Schl?mann, Journal of Applied Physics, 1962


% n_helix: Helix index (e.g. 1-10)
% s_helix : Coordinates of the point inside the rod in cylinder coordinates
% s_helix = [phi_i, r_i, theta_i];


%
clear all
close all
clc
n_helix = 6;
s_helix = [pi, 0, 0];


[h_helix, D_helix, k_helix, d_helix] = Helixinfo(n_helix);

R_helix = D_helix/2;
r_helix = d_helix/2;

phi_i = s_helix(1);
r_i = s_helix(2);
theta_i = s_helix(3);

syms s_str s phi_str r_str theta_str phi r theta

s_str = [cos(phi_str) -sin(phi_str) 0;...
    sin(phi_str) cos(phi_str) 0;...
    0 0 1]*[R_helix + r_str*cos(theta_str); 0; h_helix*phi_str/(2*pi) + r_str*sin(theta_str)];
s = [cos(phi) -sin(phi) 0;...
    sin(phi) cos(phi) 0;...
    0 0 1]*[R_helix + r*cos(theta); 0; h_helix*phi/(2*pi) + r*sin(theta)]
%%
fun_mantel = (s_str-s)*cross(diff(s_str,phi_str),diff(s_str,theta_str))'/norm(s_str-s)^3;
fun_deckel = (s_str-s)*cross(diff(s_str,phi_str),diff(s_str,r_str))'/norm(s_str-s)^3;


phi_min = 0;
phi_max = k_helix*2*pi;
theta_min = 0;
theta_max = 2*pi;
r_min = 0;
r_max = r_helix;



N_matrix = zeros(3,3);

for i = 1 : 3
    
    for k = 1 : 3
        
        
        fun_mantel_param(phi,phi_str,r,r_str,theta,theta_str) = fun_mantel(i,k);
        fun_deckel_param(phi,phi_str,r,r_str,theta,theta_str) = fun_deckel(i,k);
        
        N_mantel = 0;
        N_deckel_up = 0;
        N_deckel_down = 0;
        integrand_deckel_up = matlabFunction(fun_deckel_param(phi_i,phi_max,r_i,r_str,theta_i,theta_str));
        integrand_deckel_down = matlabFunction(fun_deckel_param(phi_i,phi_min,r_i,r_str,theta_i,theta_str));
        integrand_mantel = matlabFunction(fun_mantel_param(phi_i,phi_str,r_i,r_helix,theta_i,theta_str));
        
        
        
        if(strcmp(func2str(integrand_mantel), '@()0.0') )
            disp('up & down')
            N_deckel_up = 1/(4*pi)*abs(integral2(integrand_deckel_up,phi_min,phi_max,r_min,r_max));
            N_deckel_down = 1/(4*pi)*abs(integral2(integrand_deckel_down,phi_min,phi_max,r_min,r_max));
            
            
        elseif(strcmp(func2str(integrand_deckel_up), '@()0.0') && strcmp(func2str(integrand_deckel_down), '@()0.0'))
            disp('mantel')
            N_mantel = 1/(4*pi)*abs(integral2(integrand_mantel,phi_min,phi_max,theta_min,theta_max));
            
        elseif(strcmp(func2str(integrand_deckel_up), '@()0.0'))
            disp('mantel & down')
            N_mantel = 1/(4*pi)*abs(integral2(integrand_mantel,phi_min,phi_max,theta_min,theta_max));
            N_deckel_down = 1/(4*pi)*abs(integral2(integrand_deckel_down,phi_min,phi_max,r_min,r_max));
       
        elseif(strcmp(func2str(integrand_deckel_down), '@()0.0'))
            disp('mantel & up')
            N_mantel = 1/(4*pi)*abs(integral2(integrand_mantel,phi_min,phi_max,theta_min,theta_max));
            N_deckel_up = 1/(4*pi)*abs(integral2(integrand_deckel_up,phi_min,phi_max,r_min,r_max));
        
        else
            disp('all')
            N_mantel = 1/(4*pi)*abs(integral2(integrand_mantel,phi_min,phi_max,theta_min,theta_max));
            N_deckel_up = 1/(4*pi)*abs(integral2(integrand_deckel_up,phi_min,phi_max,r_min,r_max));
            N_deckel_down = 1/(4*pi)*abs(integral2(integrand_deckel_down,phi_min,phi_max,r_min,r_max));
            
        end
        
        N = (N_mantel + N_deckel_up + N_deckel_down);
        
        N_matrix(i,k) = N;
        
    end
end

trace = N_matrix(1,1) + N_matrix(2,2) + N_matrix(3,3)




