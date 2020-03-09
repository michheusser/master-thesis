function N_matrix = DemagFactor_Rod_Analytical(L_rod,R_rod,s_rod)
% Calculates the Demagnetization Matrix at a specific point inside a rod of
% length L_rod and radius R_rod centered at the origin and aligned with the
% z-axis.

%Source: "A Sum Rule Concerning the Inhomogeneous Demagnetizing Field in
%Nonellipsoidal Samples, Ernst Schl?mann, Journal of Applied Physics, 1962


% L : Length of the rod
% R : Radius of the rod
% s : Coordinates of the point inside the rod in cylinder coordinates
% s = [r_i, phi_i, z_i];

r_i = s_rod(1);
phi_i = s_rod(2);
z_i = s_rod(3);


syms s_str s phi_str r_str z_str phi r z

s_str = [r_str*cos(phi_str); r_str*sin(phi_str); z_str];
s = [r*cos(phi); r*sin(phi); z];


fun_mantel = (s_str-s)*cross(diff(s_str,phi_str),diff(s_str,z_str)).'/norm(s_str-s)^3;
fun_deckel = (s_str-s)*cross(diff(s_str,phi_str),diff(s_str,r_str)).'/norm(s_str-s)^3;

z_min = -L_rod/2;
z_max = L_rod/2;
phi_min = 0;
phi_max = 2*pi;
r_min = 0;
r_max = R_rod;



N_matrix = zeros(3,3);

for i = 1 : 3
    
    for k = 1 : 3
        
        
        fun_mantel_param(phi,phi_str,r,r_str,z,z_str) = fun_mantel(i,k);
        fun_deckel_param(phi,phi_str,r,r_str,z,z_str) = fun_deckel(i,k);
        
        N_mantel = 0;
        N_deckel_up = 0;
        N_deckel_down = 0;
        
        if(fun_mantel_param == 0 )
            %disp('deckel')
            integrand_deckel_up = matlabFunction(fun_deckel_param(phi_i,phi_str,r_i,r_str,z_i,L_rod/2));
            integrand_deckel_down = matlabFunction(fun_deckel_param(phi_i,phi_str,r_i,r_str,z_i,-L_rod/2));
            N_deckel_up = 1/(4*pi)*abs(integral2(integrand_deckel_up,phi_min,phi_max,r_min,r_max));
            N_deckel_down = 1/(4*pi)*abs(integral2(integrand_deckel_down,phi_min,phi_max,r_min,r_max));
            
        elseif(fun_deckel_param == 0)
            %disp('mantel')
            integrand_mantel = matlabFunction(fun_mantel_param(phi_i,phi_str,r_i,R_rod,z_i,z_str));
            N_mantel = 1/(4*pi)*abs(integral2(integrand_mantel,phi_min,phi_max,z_min,z_max));
            
        else
            disp('both')
            integrand_deckel_up = matlabFunction(fun_deckel_param(phi_i,phi_str,r_i,r_str,z_i,L_rod/2));
            integrand_deckel_down = matlabFunction(fun_deckel_param(phi_i,phi_str,r_i,r_str,z_i,-L_rod/2));
            integrand_mantel = matlabFunction(fun_mantel_param(phi_i,phi_str,r_i,R_rod,z_i,z_str));
            N_mantel = 1/(4*pi)*abs(integral2(integrand_mantel,phi_min,phi_max,z_min,z_max));
            N_deckel_up = 1/(4*pi)*abs(integral2(integrand_deckel_up,phi_min,phi_max,r_min,r_max));
            N_deckel_down = 1/(4*pi)*abs(integral2(integrand_deckel_down,phi_min,phi_max,r_min,r_max));
            
        end
        
        N = (N_mantel + N_deckel_up + N_deckel_down);
        
        N_matrix(i,k) = N;
    
    end
end

trace = N_matrix(1,1) + N_matrix(2,2) + N_matrix(3,3)

end



