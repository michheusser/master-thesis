function H_demag = DemagField_Wire_Analytical(coord, shape, shape_param, M_fun_x, M_fun_y,M_fun_z)

% Calculates de Demagnetization Field H_d at a certain point of a magnetic
% shape given the magnetization inside the body.
% M is given as

H_demag_x = 0;
H_demag_y = 0;
H_demag_z = 0;

syms s s_str s_param phi phi_str phi_param R R_str R_param

M_fun_xyz = [sym(M_fun_x);sym(M_fun_y);sym(M_fun_z)];

phi_i = coord(1);
R_i = coord(2);
s_i = coord(3);


if(strcmp(shape,'rod'))
    
    R_wire = shape_param(1);
    L_rod = shape_param(2);
    f_fun = [0; 0; s_param];
    dfds = diff(f_fun,s_param);
    z_w = dfds/norm(dfds);
    x_w = [1 0 0]';
    y_w = [0 1 0]';
    rot_Matrix(s_param) = [x_w y_w z_w];
    rot_Matrix_ = [x_w y_w z_w];
    
    phi_min = 0;
    phi_max = 2*pi;
    R_min = R_wire*0.01;
    R_max = R_wire*0.99;
    s_min = L_rod*0.01;
    s_max = L_rod*0.99;

end

assume(0 <= phi_param < 2*pi);
assume(0 <= phi_str < 2*pi);
assume(0 <= phi < 2*pi);
assume(R_min <= R_param <= R_wire);
assume(R_min <= R_str <= R_wire);
assume(R_min <= R <= R_wire);
assume(s_min <= s_param  <= s_max);
assume(s_min <= s_str <= s_max);
assume(s_min <= s <= s_max);

f_fun_x(s_param) = f_fun(1);
f_fun_y(s_param) = f_fun(2);
f_fun_z(s_param) = f_fun(3);

% M_fun_x(phi_param,R_param,s_param) = M_fun_xyz(1);
% M_fun_y(phi_param,R_param,s_param) = M_fun_xyz(2);
% M_fun_z(phi_param,R_param,s_param) = M_fun_xyz(3);

M_fun_phiRs = rot_Matrix_'*M_fun_xyz;

M_fun_phi(phi_param,R_param,s_param) = M_fun_phiRs(1);
M_fun_R(phi_param,R_param,s_param) = M_fun_phiRs(2);
M_fun_s(phi_param,R_param,s_param) = M_fun_phiRs(3);

r_str = rot_Matrix(s_str)*[R_str*cos(phi_str); R_str*sin(phi_str); 0]+[f_fun_x(s_str); f_fun_y(s_str); f_fun_z(s_str)];
r = rot_Matrix(s)*[R*cos(phi); R*sin(phi); 0]+[f_fun_x(s); f_fun_y(s);f_fun_z(s)];


dr_dR = diff(r_str,R_str);
dr_dphi = diff(r_str,phi_str);
dr_ds = diff(r_str, s_str);
h_R = simplify(sqrt(dr_dR(1)^2 + dr_dR(2)^2 + dr_dR(3)^2));
h_phi = simplify(sqrt(dr_dphi(1)^2 + dr_dphi(2)^2 + dr_dphi(3)^2));
h_s = simplify(sqrt(dr_ds(1)^2 + dr_ds(2)^2 + dr_ds(3)^2));

divM = simplify(1/(h_R*h_phi*h_s)*(diff(M_fun_R(phi_str,R_str,s_str)*h_phi*h_s,R_str) + diff(M_fun_phi(phi_str,R_str,s_str)*h_R*h_s,phi_str) + diff(M_fun_s(phi_str,R_str,s_str)*h_phi*h_R,s_str)))

jacobidet = abs(simplify(det(jacobian(r_str,[phi_str, R_str, s_str]))));

integrand_vec = (r_str-r)./((r_str(1)-r(1))^2+(r_str(2)-r(2))^2+(r_str(3)-r(3))^2)^(3/2)*divM*jacobidet;

integrand_x_param(phi, phi_str, R, R_str, s, s_str) = integrand_vec(1);
integrand_y_param(phi, phi_str, R, R_str, s, s_str) = integrand_vec(2);
integrand_z_param(phi, phi_str, R, R_str, s, s_str) = integrand_vec(3);

integrand_x(phi_str,R_str,s_str) = simplify(integrand_x_param(phi_i,phi_str,R_i,R_str,s_i,s_str));
integrand_y(phi_str,R_str,s_str) = simplify(integrand_y_param(phi_i,phi_str,R_i,R_str,s_i,s_str));
integrand_z(phi_str,R_str,s_str) = simplify(integrand_z_param(phi_i,phi_str,R_i,R_str,s_i,s_str));

integration_Method = 'iterated';

if(integrand_x(phi_str,R_str,s_str) ~= 0)
integrand_x_handle = matlabFunction(integrand_x(phi_str, R_str, s_str),'vars',[s_str R_str phi_str]);
H_demag_x = integral3(integrand_x_handle, s_min,s_max,R_min,R_max,phi_min,phi_max,'Method',integration_Method);
end

if(integrand_y(phi_str,R_str,s_str) ~= 0)
integrand_y_handle = matlabFunction(integrand_y(phi_str, R_str, s_str),'vars',[s_str R_str phi_str]);
H_demag_y = integral3(integrand_y_handle, s_min,s_max,R_min,R_max,phi_min,phi_max,'Method',integration_Method);
end

if(integrand_z(phi_str,R_str,s_str) ~= 0)
integrand_z_handle = matlabFunction(integrand_z(phi_str, R_str, s_str),'vars',[s_str R_str phi_str]);
H_demag_z = integral3(integrand_z_handle, s_min,s_max,R_min,R_max,phi_min,phi_max,'Method',integration_Method);
end


H_demag = [H_demag_x; H_demag_y; H_demag_z];

end
