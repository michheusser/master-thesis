close all
clear all
clc

syms R phi s R_h h_h R_w p s_str h_h2pi
assume(h_h2pi >=0 )
assume(R >= 0)
assume(phi >= 0)
assume(s>= 0)
assume(s_str >= 0)
assume(R_h>= 0)
assume(h_h>= 0)
assume(p>= 0)
assume(R_w >= 0)



n_helix = 3;
[h_helix, D_helix, k_helix, d_helix] = Helixinfo(n_helix);

r_c = [R_h*cos(s); R_h*sin(s); s*h_h2pi];
r_c_str = [R_h*cos(s_str); R_h*sin(s_str); s_str*h_h2pi];
x_w = [cos(s); sin(s); 0];
dr_c_ds = diff(r_c,s);
dr_c_str_ds = diff(r_c_str,s_str);
z_w = simplify(dr_c_ds/sqrt(dr_c_ds(1)^2 + dr_c_ds(2)^2 + dr_c_ds(3)^2));
y_w = cross(z_w,x_w);

A = [x_w y_w z_w];
r = A*[R*cos(phi); R*sin(phi); 0] + r_c;
%%

J = simplify(jacobian(r,[R, phi, s]));
detJ = simplify(det(J));

Integral = simplify(int(simplify(int(detJ,phi,0, 2*p)),R,0,R_w));
Integral_param(R_w,h_h2pi,R_h,p) = Integral;
Result(R_w,h_h2pi,R_h,p) = (1/(4*p))*Integral_param(R_w,h_h2pi,R_h,p);
Result_spec = double(Result(d_helix/2, h_helix/(2*pi), D_helix/2,pi));

rcstr_rc = r_c_str - r_c;
absrcstr_rc = sqrt(rcstr_rc(1)^2 + rcstr_rc(2)^2 + rcstr_rc(3)^2);
f_approx = simplify(1/(absrcstr_rc)^3*(3/absrcstr_rc^2*rcstr_rc*rcstr_rc.'-eye(3)));

syms a_0 a_1 a_2 a_3 a_4 b_0 b_1 b_2 b_3 b_4 c_0 c_1 c_2 c_3 c_4

M_s = [a_0 + a_1*s + a_2*s^2 + a_3*s^3 + a_4*s^4; b_0 + b_1*s + b_2*s^2 + b_3*s^3 + b_4*s^4; c_0 + c_1*s + c_2*s^2 + c_3*s^3 + c_4*s^4];
M_s_str = [a_0 + a_1*s_str + a_2*s_str^2 + a_3*s_str^3 + a_4*s_str^4; b_0 + b_1*s_str + b_2*s_str^2 + b_3*s_str^3 + b_4*s_str^4; c_0 + c_1*s_str + c_2*s_str^2 + c_3*s_str^3 + c_4*s_str^4];

N_op = Result_spec*f_approx*M_s_str;

%%
cross_mantel_param(R,phi,s) = simplify(cross(diff(r,phi),diff(r,s)))
cross_mantel = cross_mantel_param(R_w,phi,s)
