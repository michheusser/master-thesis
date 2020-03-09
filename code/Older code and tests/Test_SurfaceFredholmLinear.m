close all
clear all
clc

n_helix = 3;
[h_helix, D_helix, k_helix, d_helix] = Helixinfo(n_helix);

s_min = 0;
s_max = 2*pi*k_helix;

L = 1;
N = 5;


syms x y t
U_exact_sym_param(t) = [2*t; 3*t^2; -t^3];
%K_sym =  [cos(y)*x x^2*y x*exp(y/x); x*y^2 exp(x*y)/y x+y*log(x*y); cos(x)*cos(y) sin(y)*exp(cos(x)*y) x^2*cos(y)];
K_sym = [x*y^2 x^2*y y^3; x^2*y x*y x*y^2; y^3 x*y^2 x^2*y];

U_exact_sym = U_exact_sym_param(x);
U_exact_sym_y = U_exact_sym_param(y);

F_sym =  U_exact_sym - int(K_sym*U_exact_sym_y,y,a,b);

K_handle = matlabFunction(K_sym);
F_handle = matlabFunction(F_sym);

U_approx_handle = Fredholm_SurfaceLinear(K_handle,F_handle,L,a,b,N,c);


U_approx_sym = sym(U_approx_handle)

close all
s = length(U_exact_sym);
for i = 1 : s
    subplot(s,1,i)
    ezplot(U_exact_sym(i),[a,b])
    hold on
    ezplot(U_approx_sym(i),[a,b])
    legend(['U_{' num2str(i) ',exact}(x)'],['U_{' num2str(i) ',approx}(x)'])
    hold off
end




