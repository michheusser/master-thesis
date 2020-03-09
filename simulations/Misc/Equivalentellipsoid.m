close all
clear all
clc

syms x y
P = vpasolve([x^3 + 2*x == y, y^2 == x], [x, y])

P.x
P.y
%%


% m = 0.5
% n = 0.3
% (m*n)/(sqrt(1-n^2)*(1-m^2))*(ellipticF(acos(n),(1-m^2)/(1-n^2))-ellipticE(acos(n),(1-m^2)/(1-n^2)))
% (m*n)/((m^2-n^2)*sqrt(1-n^2))*(((sqrt(1-n^2)*m)/n)-(ellipticE(acos(n),(1-m^2)/(1-n^2))))

N_11 = 0.33
N_22 = 0.33

syms m n
S = vpasolve([(m*n)/(sqrt(1-n^2)*(1-m^2))*(ellipticF(acos(n),(1-m^2)/(1-n^2))-ellipticE(acos(n),(1-m^2)/(1-n^2))) == N_11, (m*n)/((m^2-n^2)*sqrt(1-n^2))*(((sqrt(1-n^2)*m)/n)-(ellipticE(acos(n),(1-m^2)/(1-n^2)))) == N_22], [m, n])
%%
S.m
S.n