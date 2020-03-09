%function f = Fredholm_Homotopymethod(K_handle,g_handle,a,b,N)
%%
close all
clear all
clc

syms x t p f

f_t = [t; t^2];
f_x = [x; x^2];
K = [x^2 t*x; t*x t^2];
g = f_x - int(K*f_t,t,0,1);

g_handle = matlabFunction(g);
K_handle = matlabFunction(K);

%K_handle = @(x,t) [x^2 t*x; t*x t^2];
%g_handle = @(x) [ x - (x*(2*x + 1))/4; x^2 - x/3 - 1/5];

a = 0;
b = 1;
N = 20;

% Solves the equation:
% f(x) = g(x) + int_a^b K(x,t,f(t)) dt
% Where f and g are a column vector of dimension dim and K a
% nonlinear operator

% Inputs have to be function handles:


K = sym(K_handle);
g = sym(g_handle);
assert(length(K)==length(g))
dim = length(g);

p_vec = zeros(N+1,1)*p;
for i = 0 : N
    
    p_vec(i+1,1) = p^(i);
    
end

f = zeros(dim,N+1)*x;
h = zeros(dim,N)*x;

n = 0;
f(:,n+1) = g;
h(:,n+1) = f(:,n+1);
n = 1;
f(:,n+1) = -int(K*h(:,(n+1)-1),t,a,b);

f_sol = f(:,0+1) + f(:,1+1);
for n = 2 : N
    
    h(:,(n+1)-1) = f*p_vec;
    f(:,(n+1)) = -int(K*(h(:,(n+1)-1)-h(:,(n+1)-2)),t,a,b);
    f_sol = f_sol + f(:,(n+1));
end

f_sol_param(p,x) = simplify(f_sol);
f_sol_final = simplify(limit(f_sol_param,p,1));

f_sol_final_1 = simplify(limit(f_sol(1),p,1));
f_sol_final_2 = simplify(limit(f_sol(2),p,1));


ezplot(f_sol_final_1,[a,b])
hold on
ezplot(f_sol_final_2,[a,b])
hold off