%function U_sol_handle = Fredholm_Linear_Homotopy(K_handle,F_handle,L,a,b,N)

% Solves the Fredholm Linear System of Equations of the Second Kind:
% U(x) = F(x) + lambda*Integral_a^b(K(x,y)*U(y))
%   
% Where U(x) is the solution to be calculated
%
% K : K(x,y) Kernel
% F : F(x) Function
% L : Lambda parameter
% [a,b]: Limits of the integral
% N: Number of taylor terms
% c: Point in [a,b] where the taylor expansion is done
%
% Source: "Numerical Solution of the System of Linear Fredholm Integral
% Equations Based on Degenerating Kernels"
% S. Karimi, M. Jozi
% TWMS J. Pure Appl. Math V.6, N.1, 2015, pp. 111-119
%
% Author: Michel Heusser (ETH Z?rich)
%

close all
clear all
clc

a = 3;
b = 5;
c = (a+b)/2;
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


syms x y

K = sym(K_handle);
F = sym(F_handle);
assert(length(K)==length(F))
assert(c >= a && c <=b)
n = length(F);

f = cell(n,N);

f_0 = F;

for i = 1 : n
     
    for j = 1 : n
       
        for k = 0 : q(i,j)
            
        end
        
    end
    
end

