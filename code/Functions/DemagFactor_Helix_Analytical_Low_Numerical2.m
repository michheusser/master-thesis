function [N_global_cell, N_local_cell, N_average_1, N_average_2] = DemagFactor_Helix_Analytical_Low(n_helix,chi_m,N,shape_function, integration_method, s_vec)

%% Symbolic Variables
syms s s_pr phi phi_pr R R_pr
assume(s >= 0)
assume(s_pr >= 0)
assume(1 > sin(phi) >= 0)
assume(1 > sin(phi_pr) >= 0)
assume(1 > cos(phi) >= -1)
assume(1 > cos(phi_pr) >= -1)
assume(1 > sin(s) >= 0)
assume(1 > sin(s_pr) >= 0)
assume(1 > cos(s) >= -1)
assume(1 > cos(s_pr) >= -1)
assume(R >= 0)
assume(R_pr >= 0)

%% Further Variable Definitions
mu_0 = 1.25663706e-6;

[h_helix, D_helix, k_helix, d_helix] = Helixinfo(n_helix);

L = -chi_m/(4*pi);
F_sym(s) = chi_m*eye(3)+0*s;

s_min = 0;
s_max = 2*pi*k_helix;

phi_min = 0;
phi_max = 2*pi;

R_min = 0;
R_max = d_helix/2;



%% Construction of r(R,phi,s)

r_c(s) = [(D_helix/2)*cos(s); (D_helix/2)*sin(s); h_helix/(2*pi)*s];

dfds(s) = diff(r_c(s),s);
z_w(s) = dfds/norm(dfds(s));
x_w(s) = [cos(s) -sin(s) 0; sin(s) cos(s) 0; 0 0 1]*[1; 0; 0];
y_w(s) = cross(z_w(s),x_w(s));
rot_Matrix(s) = [x_w(s) y_w(s) z_w(s)];

r(R,phi,s) = rot_Matrix(s)*[R*cos(phi); R*sin(phi); 0]+r_c(s);

%% Construction of surface integral kernels
disp('Constructing general surface integral kernels')
K_tube(R,phi,s,R_pr,phi_pr,s_pr) = simplify(cross(diff(r(R_pr,phi_pr,s_pr),phi_pr),diff(r(R_pr,phi_pr,s_pr),s_pr))*(r(R_pr,phi_pr,s_pr)-r(R,phi,s)).'/norm(r(R_pr,phi_pr,s_pr)-r(R,phi,s))^3);
K_top(R,phi,s,R_pr,phi_pr,s_pr) = simplify(cross(diff(r(R_pr,phi_pr,s_pr),R_pr),diff(r(R_pr,phi_pr,s_pr),phi_pr))*(r(R_pr,phi_pr,s_pr)-r(R,phi,s)).'/norm(r(R_pr,phi_pr,s_pr)-r(R,phi,s))^3);
K_bottom(R,phi,s,R_pr,phi_pr,s_pr) = simplify(cross(diff(r(R_pr,phi_pr,s_pr),phi_pr),diff(r(R_pr,phi_pr,s_pr),R_pr))*(r(R_pr,phi_pr,s_pr)-r(R,phi,s)).'/norm(r(R_pr,phi_pr,s_pr)-r(R,phi,s))^3);
disp('Kernels completed')



%% Construction of matrix function
disp('Constructing Basic Functions')
if(strcmp(shape_function,'legendre'))
    f_shape = zeros(N,1)*s;
    
    for i = 1 : N
        
        n= i-1;
        shift_x = (s_max-s_min)/2;
        scale_x = (s_max-s_min)/2;
        f_base(s) = 1/(2^n*factorial(n))*diff((s^2-1)^n,n);
        f_shape(i) = f_base(s);
        f_shape(i) = f_base((s-shift_x)/scale_x);
        
    end
    
    
elseif(strcmp(shape_function,'fourier'))
    N = 2*N+1;
    f_shape = zeros(N,1)*s;
    
    for n = 1 : (N-1)/2
        f_shape(2*n-1) = sin(2*pi*n*s/(2*pi*k_helix));
        f_shape(2*n) = cos(2*pi*n*s/(2*pi*k_helix));
    end
    f_shape(N) = 1;
    
end

M_pack(s) = [f_shape' zeros(1,N) zeros(1,N);...
    zeros(1,N) f_shape' zeros(1,N);...
    zeros(1,N) zeros(1,N) f_shape'];
disp('Basic Functions completed')

%% Construction of integrand function handles

disp('Constructing function handles')
integrand_tube = K_tube(R,phi,s,R_pr,phi_pr,s_pr)*M_pack(s_pr);
integrand_top = K_top(R,phi,s,R_pr,phi_pr,s_pr)*M_pack(s_pr);
integrand_bottom = K_bottom(R,phi,s,R_pr,phi_pr,s_pr)*M_pack(s_pr);
disp('Function handles completed')

[row_zero_tube, col_zero_tube] = find(integrand_tube);
[row_zero_top, col_zero_top] = find(integrand_top);
[row_zero_bottom, col_zero_bottom] = find(integrand_bottom);

disp('Constructing handles cells')
integrand_tube_cell = cell(3,3*N);
integrand_top_cell = cell(3,3*N);
integrand_bottom_cell = cell(3,3*N);

disp(['Progress: 0%'])
for i = 1 : 3
    
    for j = 1 : N*3
        
        integrand_tube_cell{i,j} = matlabFunction(integrand_tube(i,j),'vars',[R,phi,s,R_pr,phi_pr,s_pr]);
        integrand_top_cell{i,j} = matlabFunction(integrand_top(i,j),'vars',[R,phi,s,R_pr,phi_pr,s_pr]);
        integrand_bottom_cell{i,j} = matlabFunction(integrand_bottom(i,j),'vars',[R,phi,s,R_pr,phi_pr,s_pr]);
        disp(['Progress: ' num2str(round( 100*((i-1)*N*3+j)/(3*3*N) )) '%'])
    end
end
disp('Cells handles completed')

%% Construction of Integrands and Integrals
R_i = 0;
phi_i = 0;
s_i_vec = 2*pi*k_helix*(1:N)'/(N+1);

A_system = zeros(3*N);

disp(['Calculating Approximation of order N = ' num2str(N) ' using ' shape_function])
for j = 1 : N
    
    disp(['Starting Step: ' num2str(j) '/' num2str(N)])
    s_i = s_i_vec(j);
    Integral_Matrix_Tube = zeros(3,3*N);
    Integral_Matrix_Top = zeros(3,3*N);
    Integral_Matrix_Bottom = zeros(3,3*N);
    
    for i = 1 : length(row_zero_tube)
        m = row_zero_tube(i);
        n = col_zero_tube(i);
        integrand_tube_handle = integrand_tube_cell{m,n};
        integrand_tube_handle_specific = @(phi_pr,s_pr) integrand_tube_handle(R_i,phi_i,s_i,R_max,phi_pr,s_pr);
        integral_tube = integral2(integrand_tube_handle_specific,phi_min,phi_max,s_min,s_max,'Method',integration_method);
        Integral_Matrix_Tube(m,n) = integral_tube;
        disp(['Progress: Tube Integrals ' num2str(round( 100*(i/length(row_zero_tube)))) '%'])
    end
    
    for i = 1 : length(row_zero_top)
        m = row_zero_top(i);
        n = col_zero_top(i);
        integrand_top_handle = integrand_tube_cell{m,n};
        integrand_top_handle_specific = @(R_pr,phi_pr) integrand_top_handle(R_i,phi_i,s_i,R_pr,phi_pr,s_max);
        integral_top = integral2(integrand_top_handle_specific,R_min,R_max,phi_min,phi_max,'Method',integration_method);
        Integral_Matrix_Top(m,n) = integral_top;
        disp(['Progress: Top Integrals ' num2str(round( 100*(i/length(row_zero_top)))) '%'])
    end
    
    for i = 1 : length(row_zero_bottom)
        m = row_zero_bottom(i);
        n = col_zero_bottom(i);
        integrand_bottom_handle = integrand_tube_cell{m,n};
        integrand_bottom_handle_specific = @(R_pr,phi_pr) integrand_bottom_handle(R_i,phi_i,s_i,R_pr,phi_pr,s_max);
        integral_bottom = integral2(integrand_bottom_handle_specific,R_min,R_max,phi_min,phi_max,'Method',integration_method);
        Integral_Matrix_Bottom(m,n) = integral_bottom;
        disp(['Progress: Bottom Integrals ' num2str(round( 100*(i/length(row_zero_bottom)))) '%'])
    end
    
    
    Integral_Matrix = chi_m/(4*pi)*(Integral_Matrix_Tube + Integral_Matrix_Top + Integral_Matrix_Bottom);
    
    A_system(3*(j-1)+1:3*(j-1)+3,1:3*N) = M_pack(s_i) + Integral_Matrix;
    b_system(3*(j-1)+1:3*(j-1)+3,1:3) = chi_m*eye(3);
    
    
    
end

disp(['Step ' num2str(j) '/' num2str(N) ' completed.'])


%% Coefficient Calculation
assignin('base', 'A_system', A_system)

PHI = A_system\b_system;
Psi_sol_global(s) = M_pack(s)*PHI;


%% N-Matrix Calculation

N_global_cell = cell(length(s_vec),1);
N_local_cell = cell(length(s_vec),1);

N_average_1 = zeros(3,3);
Psi_inv_av = zeros(3,3);
Psi_av = zeros(3,3);

for i = 1 : length(s_vec)
    
N_global_cell{i,1} = inv(double(Psi_sol_global(s_vec(i))))-eye(3)/chi_m;
N_local_cell{i,1} = double(rot_Matrix(s_vec(i)).'*N_global_cell{i,1}*rot_Matrix(s_vec(i)));

N_average_1 = N_average_1 + 1/length(s_vec)*N_global_cell{i,1};

Psi_av = Psi_av + 1/length(s_vec)*double(Psi_sol_global(s_vec(i)));
Psi_inv_av = Psi_inv_av + 1/length(s_vec)*inv(double(Psi_sol_global(s_vec(i))));
end

N_average_2 = Psi_inv_av - 1/chi_m*Psi_inv_av*Psi_av;


 