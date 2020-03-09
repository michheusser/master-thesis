function [N_matrix_global,N_matrix_local, rot_Matrix_local] = DemagFactor_Wire_Analytical(coord,shape,shape_param)

% Calculates the Demagnetization Matrix at a specific point inside a helix
% of helix index n_helix

%Source: "A Sum Rule Concerning the Inhomogeneous Demagnetizing Field in
%Nonellipsoidal Samples, Ernst Schl?mann, Journal of Applied Physics, 1962


% shape: 'helix', 'rod', 'other'
% Helix: shape_param = n_helix
% Thin Helix: shape_param = n_helix
% Rod: shape_param = [R_wire, L_rod]
% Other: shape_param = [f_fun(s_param), R_wire, s_min, s_max, R_rot(s_param)]

syms s_param s s_str phi phi_str R R_str

phi_i = coord(1);
R_i = coord(2);
s_i = coord(3);

if(strcmp(shape,'helix'))
    
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
    
    R_min = 0;
    R_max = R_wire;
    
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
    
    R_min = 0;
    R_max = R_wire;
end


if(s_i == s_min)
    s_i = s_min+sqrt(eps);
elseif(s_i == s_max)
    s_i = s_max-sqrt(eps);
end

if(R_i == R_max)
    R_i = R_max-sqrt(eps);
elseif(R_i == R_min)
    R_i = R_min + sqrt(eps);
end


f_fun_x(s_param) = f_fun(1);
f_fun_y(s_param) = f_fun(2);
f_fun_z(s_param) = f_fun(3);

r_str = rot_Matrix(s_str)*[R_str*cos(phi_str); R_str*sin(phi_str); 0]+[f_fun_x(s_str); f_fun_y(s_str); f_fun_z(s_str)];
r = rot_Matrix(s)*[R*cos(phi); R*sin(phi); 0]+[f_fun_x(s); f_fun_y(s);f_fun_z(s)];

fun_mantel = (r_str-r)*cross(diff(r_str,phi_str),diff(r_str,s_str)).'/((r_str(1)-r(1))^2+(r_str(2)-r(2))^2+(r_str(3)-r(3))^2)^(3/2);
fun_deckel_up = (r_str-r)*cross(diff(r_str,R_str),diff(r_str,phi_str)).'/((r_str(1)-r(1))^2+(r_str(2)-r(2))^2+(r_str(3)-r(3))^2)^(3/2);
fun_deckel_down = (r_str-r)*cross(diff(r_str,phi_str),diff(r_str,R_str)).'/((r_str(1)-r(1))^2+(r_str(2)-r(2))^2+(r_str(3)-r(3))^2)^(3/2);
fun_mantel_in = (r_str-r)*cross(diff(r_str,s_str),diff(r_str,phi_str)).'/((r_str(1)-r(1))^2+(r_str(2)-r(2))^2+(r_str(3)-r(3))^2)^(3/2);

if(~strcmp(shape, 'thin helix'))
    fun_mantel_in = fun_mantel_in*zeros(3,3);
end


phi_min = 0;
phi_max = 2*pi;


integration_method = 'iterated';

N_matrix_global = zeros(3,3);
for i = 1 : 3
    
    for k = 1 : 3
        
        disp(['Calculating Tensor Component N_' num2str(i) num2str(k)])
        fun_mantel_param(phi,phi_str,R,R_str,s,s_str) = fun_mantel(i,k);
        fun_mantel_in_param(phi,phi_str,R,R_str,s,s_str) = fun_mantel_in(i,k);
        fun_deckel_up_param(phi,phi_str,R,R_str,s,s_str) = fun_deckel_up(i,k);
        fun_deckel_down_param(phi,phi_str,R,R_str,s,s_str) = fun_deckel_down(i,k);
        
        %simplify(fun_deckel_down_param(phi,phi_str,R,R_str,s,s_str))
        
        integrand_deckel_up_param(phi_str,R_str)= simplify(fun_deckel_up_param(phi_i,phi_str,R_i,R_str,s_i,s_max));
        integrand_deckel_down_param(phi_str,R_str)= simplify(fun_deckel_down_param(phi_i,phi_str,R_i,R_str,s_i,s_min));
        integrand_mantel_param(phi_str,s_str) = simplify(fun_mantel_param(phi_i,phi_str,R_i,R_max,s_i,s_str));
        integrand_mantel_in_param(phi_str,s_str) = simplify(fun_mantel_in_param(phi_i,phi_str,R_i,R_min,s_i,s_str));
        
        N_mantel = 0;
        N_mantel_in = 0;
        N_deckel_up = 0;
        N_deckel_down = 0;
        
        if(integrand_deckel_up_param(phi_str,R_str) ~= 0)
            %disp('up')
            integrand_deckel_up_handle = matlabFunction(integrand_deckel_up_param(phi_str,R_str),'vars', [R_str phi_str]);
            N_deckel_up = 1/(4*pi)*(integral2(integrand_deckel_up_handle,R_min,R_max,phi_min,phi_max,'Method',integration_method));
        end
        
        if(integrand_deckel_down_param(phi_str,R_str) ~= 0)
            %disp('down')
            integrand_deckel_down_handle = matlabFunction(integrand_deckel_down_param(phi_str,R_str),'vars', [R_str phi_str]);
            matlabFunction(integrand_deckel_down_param(phi_str,R_str),'vars', [R_str phi_str]);
            N_deckel_down = 1/(4*pi)*(integral2(integrand_deckel_down_handle,R_min,R_max,phi_min,phi_max,'Method',integration_method));
        end
        
        if(integrand_mantel_param(phi_str,s_str) ~= 0)
            %disp('mantel')
            integrand_mantel_handle = matlabFunction(integrand_mantel_param(phi_str,s_str),'vars', [phi_str s_str]);
            N_mantel = 1/(4*pi)*(integral2(integrand_mantel_handle,phi_min,phi_max,s_min,s_max,'Method',integration_method));
        end
        
        if(integrand_mantel_in_param(phi_str,s_str) ~= 0)
            %disp('mantel in')
            integrand_mantel_in_handle = matlabFunction(integrand_mantel_in_param(phi_str,s_str),'vars', [phi_str s_str]);
            N_mantel_in = 1/(4*pi)*(integral2(integrand_mantel_in_handle,phi_min,phi_max,s_min,s_max,'Method',integration_method));
        end
        
        N = (N_mantel + N_mantel_in + N_deckel_up + N_deckel_down);
        N_matrix_global(i,k) = N;
        
    end
end


rot_Matrix_local = double(rot_Matrix(s_i));
N_matrix_local = rot_Matrix_local'*N_matrix_global*rot_Matrix_local;



% %% Shows Normal Vectors
% 
% r_x(R,phi,s) = r(1);
% r_y(R,phi,s) = r(2);
% r_z(R,phi,s) = r(3);
% 
% n_normal_mantel = cross(diff(r,phi),diff(r,s))/norm(cross(diff(r,phi),diff(r,s)));
% n_normal_mantel_x(R,phi,s) = n_normal_mantel(1);
% n_normal_mantel_y(R,phi,s) = n_normal_mantel(2);
% n_normal_mantel_z(R,phi,s) = n_normal_mantel(3);
% 
% n_normal_up = cross(diff(r,R),diff(r,phi))/norm(cross(diff(r,R),diff(r,phi)));
% n_normal_up_x(R,phi,s) = n_normal_up(1);
% n_normal_up_y(R,phi,s) = n_normal_up(2);
% n_normal_up_z(R,phi,s) = n_normal_up(3);
% 
% n_normal_down = cross(diff(r,phi),diff(r,R))/norm(cross(diff(r,phi),diff(r,R)));
% n_normal_down_x(R,phi,s) = n_normal_down(1);
% n_normal_down_y(R,phi,s) = n_normal_down(2);
% n_normal_down_z(R,phi,s) = n_normal_down(3);
% 
% 
% N_s = 10;
% N_r = 10;
% N_phi = 10;
% 
% phi_disp = [0 : N_phi-1]'/N_phi*(phi_max-phi_min) + phi_min;
% R_disp = [0 : (N_r-1)]'/(N_r-1) * (R_max - R_min) + R_min;
% s_disp = [0 : (N_s-1)]'/(N_s-1) * (s_max - s_min) + s_min;
% 
% for ii = 1 : length(s_disp)
%     for jj = 1 : length(phi_disp)
%         
%         R_disp_i = R_disp(end);
%         s_disp_i = s_disp(ii);
%         phi_disp_i = phi_disp(jj);
%         
%         x_mantel(ii,jj) = r_x(R_disp_i,phi_disp_i,s_disp_i);
%         y_mantel(ii,jj) = r_y(R_disp_i,phi_disp_i,s_disp_i);
%         z_mantel(ii,jj) = r_z(R_disp_i,phi_disp_i,s_disp_i);
%         x_mantel_normal(ii,jj) = n_normal_mantel_x(R_disp_i,phi_disp_i,s_disp_i);
%         y_mantel_normal(ii,jj) = n_normal_mantel_y(R_disp_i,phi_disp_i,s_disp_i);
%         z_mantel_normal(ii,jj) = n_normal_mantel_z(R_disp_i,phi_disp_i,s_disp_i);
%         
%     
%     end
% end
% 
% for ii = 1 : length(phi_disp)
%     for jj = 1 : length(R_disp)
%         
%         R_disp_i = R_disp(jj);
%         s_disp_i = s_disp(end);
%         phi_disp_i = phi_disp(ii);
%         
%         x_up(ii,jj) = r_x(R_disp_i,phi_disp_i,s_disp_i);
%         y_up(ii,jj) = r_y(R_disp_i,phi_disp_i,s_disp_i);
%         z_up(ii,jj) = r_z(R_disp_i,phi_disp_i,s_disp_i);
%         x_up_normal(ii,jj) = n_normal_up_x(R_disp_i,phi_disp_i,s_disp_i);
%         y_up_normal(ii,jj) = n_normal_up_y(R_disp_i,phi_disp_i,s_disp_i);
%         z_up_normal(ii,jj) = n_normal_up_z(R_disp_i,phi_disp_i,s_disp_i);
%         
%     
%     end
% end
% 
% for ii = 1 : length(phi_disp)
%     for jj = 1 : length(R_disp)
%         
%         R_disp_i = R_disp(jj);
%         s_disp_i = s_disp(1);
%         phi_disp_i = phi_disp(ii);
%         
%         x_down(ii,jj) = r_x(R_disp_i,phi_disp_i,s_disp_i);
%         y_down(ii,jj) = r_y(R_disp_i,phi_disp_i,s_disp_i);
%         z_down(ii,jj) = r_z(R_disp_i,phi_disp_i,s_disp_i);
%         x_down_normal(ii,jj) = n_normal_down_x(R_disp_i,phi_disp_i,s_disp_i);
%         y_down_normal(ii,jj) = n_normal_down_y(R_disp_i,phi_disp_i,s_disp_i);
%         z_down_normal(ii,jj) = n_normal_down_z(R_disp_i,phi_disp_i,s_disp_i);
%         
%     
%     end
% end
% 
% N_fun_points = 50;
% s_points = s_min + (s_max - s_min)*[0:N_fun_points]/N_fun_points
% f_fun_x(s_param) = f_fun(1);
% f_fun_y(s_param) = f_fun(2);
% f_fun_z(s_param) = f_fun(3);
% 
% 
% x_points = zeros(length(s_points),1);
% y_points = zeros(length(s_points),1);
% z_points = zeros(length(s_points),1);
% 
% for i = 1 : length(s_points)
%     x_points(i) = double(f_fun_x(s_points(i)));
%     y_points(i) = double(f_fun_y(s_points(i)));
%     z_points(i) = double(f_fun_z(s_points(i)));
% end
% 
% plot3(x_points,y_points,z_points)
% hold on
% quiver3(x_mantel,y_mantel,z_mantel,x_mantel_normal,y_mantel_normal,z_mantel_normal,'AutoScale','on','Color','red')
% quiver3(x_up,y_up,z_up,x_up_normal,y_up_normal,z_up_normal,'AutoScale','on','Color','red')
% quiver3(x_down,y_down,z_down,x_down_normal,y_down_normal,z_down_normal,'AutoScale','on','Color','red')
% hold off
% 
% 
% 
