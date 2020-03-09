function U_sol = Fredholm_Linear(K_handle,F_handle,L,a,b,N,c)

% Source: "Numerical Solution of the System of Linear Fredholm Integral
% Equations Based on Degenerating Kernels"
% S. Karimi, M. Jozi
% TWMS J. Pure Appl. Math V.6, N.1, 2015, pp. 111-119

%close all
%clear all
%clc
%U = [x; x^2]
%K_handle = @(x,y) [x^2 y*x; y*x y^2];
%F_handle = @(x) [ x - (x*(2*x + 1))/4; x^2 - x/3 - 1/5];
%a = 0;
%b = 1;
%c = (a+b)/2;
%L = 1;
%N = 10;

syms x y
K = sym(K_handle);
F = sym(F_handle);
assert(length(K)==length(F))
s = length(F);

A_matrix_cell = cell(s,s);

      disp('0%')
for r = 1 : s
    
    for p = 1 : s
  
        A_curr = zeros(N,N);
        for i = 1 : N
            for j = 1 : N
                K_curr(x,y) = K(r,p);
                diff_K(x,y) = diff(K_curr(x,y),y,j-1);
                mu(x) = (x-c)^(i-1);
                lambda(x) = 1/factorial(j-1)*diff_K(x,c);
                A_curr(i,j) = L*int(mu(x)*lambda(x),x,a,b);
                
            end
        end
        if(r == p)
            A_matrix_cell{r,p} = eye(N) - A_curr;
        else
            A_matrix_cell{r,p} = A_curr;
        end
       disp([num2str(100*((r-1)*s+p)/(s^2)) '%'])
    end
end

F_vector_cell = cell(s,1);
for r = 1 : s
   
    F_curr = zeros(N,1);
    for j = 1 : N
        f(x) = F(r);
        mu(x) = (x-c)^(j-1);
        F_curr(j,1) = int(mu(x)*f(x),x,a,b);
    end
    
    F_vector_cell{r,1} = F_curr;
end

A_sys = zeros(N*s,N*s);
F_sys = zeros(N*s,1);

for i = 1 : s
    for j = 1 : s
        A_sys(((i-1)*N)+1:((i-1)*N)+N,((j-1)*N)+1:((j-1)*N)+N) = A_matrix_cell{i,j};
    end
    F_sys(((i-1)*N)+1:((i-1)*N)+N,1) = F_vector_cell{i,1};
end

beta = A_sys\F_sys;

A_x_matrix_cell = cell(s,s);
for r = 1 : s
   
    for p = 1 : s
       
        Lambda_curr = zeros(1,N)*x;
        for i = 1 : N
            K_curr(x,y) = K(r,p);
            diff_K(x,y) = diff(K_curr(x,y),y,i-1);
            lambda(x) = 1/factorial(i-1)*diff_K(x,c);
            Lambda_curr(1,i) = lambda(x);
        end
        A_x_matrix_cell{r,p} = Lambda_curr;
        
    end
    
end

A_x_sys = zeros(s,N)*x;
for r = 1 : s
    
   for p = 1 : s
      
       A_x_sys(r,((p-1)*N)+1:((p-1)*N)+N) = A_matrix_cell{r,p};
       
   end
    
end

U_sol = F + L*A_x_sys*beta;
U_sol_handle = matlabFunction(U_sol)