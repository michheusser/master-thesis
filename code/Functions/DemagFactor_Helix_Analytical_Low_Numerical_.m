%function N = DemagFactor_Helix_Analytical_Low(n_helix,chi_m,N,shape_function, integration_method)
N = 2;
n_helix = 3;
chi_m = 1000;
shape_function = 'fourier';
integration_method = 'auto';
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
K_tube(R,phi,s,R_pr,phi_pr,s_pr) = simplify((r(R_pr,phi_pr,s_pr)-r(R,phi,s))*cross(diff(r(R_pr,phi_pr,s_pr),phi_pr),diff(r(R_pr,phi_pr,s_pr),s_pr)).'/norm(r(R_pr,phi_pr,s_pr)-r(R,phi,s))^3);
K_top(R,phi,s,R_pr,phi_pr,s_pr) = simplify((r(R_pr,phi_pr,s_pr)-r(R,phi,s))*cross(diff(r(R_pr,phi_pr,s_pr),R_pr),diff(r(R_pr,phi_pr,s_pr),phi_pr)).'/norm(r(R_pr,phi_pr,s_pr)-r(R,phi,s))^3);
K_bottom(R,phi,s,R_pr,phi_pr,s_pr) = simplify((r(R_pr,phi_pr,s_pr)-r(R,phi,s))*cross(diff(r(R_pr,phi_pr,s_pr),phi_pr),diff(r(R_pr,phi_pr,s_pr),R_pr)).'/norm(r(R_pr,phi_pr,s_pr)-r(R,phi,s))^3);
disp('Kernels completed')
%% Construction of matrix function
if(strcmp(shape_function,'polynomials'))
    f_shape = zeros(N,1)*s;
    %for n = 1 : N
    %   f_shape(n) = (s-(s_max-s_min)/2)^(n-1);
    %end
    for n = 1 : N
        k = (n-1);
        f_shape(n) = (s-(s_max-s_min)/2)^k;
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

%% Construction of Integrands and Integrals
R_i = 0;
phi_i = 0;
s_i_vec = 2*pi*k_helix*(1:N)'/(N+1);

A_system = zeros(3*N);

integrand_tube_handle = matlabFunction(simplify(K_tube(R_i,phi_i,s,R_max,phi_pr,s_pr))*M_pack(s_pr),'vars',[phi_pr,s_pr,s]);
integrand_top_handle = matlabFunction(simplify(K_top(R_i,phi_i,s,R_pr,phi_pr,s_max))*M_pack(s_max),'vars',[R_pr,phi_pr,s]);
integrand_bottom_handle = matlabFunction(simplify(K_bottom(R_i,phi_i,s,R_pr,phi_pr,s_min))*M_pack(s_min),'vars',[R_pr,phi_pr,s]);

Integral_Matrix = zeros(3,3*N);
disp(['Calculating Approximation of order N = ' num2str(N) ' using ' shape_function])
for j = 1 : N
    disp(['Starting Step: ' num2str(j) '/' num2str(N)])
    s_i = s_i_vec(j);
    
    disp('Integrand Kernels Created')
    
    
    integrand_tube_handle_curr = @(phi_pr,s_pr) integrand_tube_handle(phi_pr,s_pr,s_i);
    integrand_top_handle_curr = @(R_pr,phi_pr) integrand_top_handle(R_pr,phi_pr,s_i);
    integrand_bottom_handle_curr = @(R_pr,phi_pr) integrand_bottom_handle(R_pr,phi_pr,s_i);
    
    integral_tube = integral(@(s_pr) integral(@(phi_pr) integrand_tube_handle_curr(phi_pr,s_pr),phi_min,phi_max,'ArrayValued',true),s_min,s_max,'ArrayValued',true);  
    integral_top = integral(@(phi_pr) integral(@(R_pr) integrand_top_handle_curr(R_pr,phi_pr),R_min,R_max,'ArrayValued',true),phi_min,phi_max,'ArrayValued',true);
    integral_bottom = integral(@(phi_pr) integral(@(R_pr) integrand_bottom_handle_curr(R_pr,phi_pr),R_min,R_max,'ArrayValued',true),phi_min,phi_max,'ArrayValued',true);
   
    %Q = integral(@(y) integral(@(x)A(x,y) ,0,2, 'ArrayValued',true), 3,5, 'ArrayValued',true)

    
%     integral_tube = integral2(@(phi_pr,s_pr) integrand_tube_handle(phi_pr,s_pr,s_i),phi_min,phi_max,s_min,s_max,'Method',integration_method);  
%     integral_top = integral2(@(R_pr,phi_pr) integrand_top_handle(R_pr,phi_pr,s_i),R_min,R_max,phi_min,phi_max,'Method',integration_method);
%     integral_bottom = integral2(@(R_pr,phi_pr) integrand_bottom_handle(R_pr,phi_pr,s_i),R_min,R_max,phi_min,phi_max,'Method',integration_method);
%     
    Integral_Matrix = chi_m/(4*pi)*(integral_tube + integral_top + integral_bottom);
      
    A_system(3*(j-1)+1:3*(j-1)+3,1:3*N) = M_pack(s_i) + Integral_Matrix;
    b_system(3*(j-1)+1:3*(j-1)+3,1:3) = chi_m*eye(3);
    
    disp(['Step ' num2str(j) '/' num2str(N) ' completed.'])
end

%% Coefficient Calculation
assignin('base', 'A_system', A_system)

PHI = A_system\b_system;

%% Magnetization Calculation

Psi_sol_global(s) = M_pack(s)*PHI;
Psi_sol_global_plot(s) = Psi_sol_global(s);
Psi_sol_global_handle = matlabFunction(Psi_sol_global(s),'vars',s);

Psi_sol_local(s) = rot_Matrix(s).'*Psi_sol_global(s)*rot_Matrix(s);
Psi_sol_local_plot = Psi_sol_local(s);
Psi_sol_local_handle = matlabFunction(Psi_sol_local(s),'vars',s);

%N_sol_av = 1/(D_helix*pi*k_helix)*integral(matlabFunction(M_sol_global(s),'vars',s),s_min,s_max,'ArrayValued',true);

%% Plot
close all
N_plot = 100;
s_plot = s_min + [1:N_plot]/(N_plot+1)*(s_max-s_min);
N_11_vec = zeros(length(s_plot),1);
N_22_vec = zeros(length(s_plot),1);
N_33_vec = zeros(length(s_plot),1);

for i = 1 : length(s_plot)
    s_i = s_plot(i);
    N_curr_global = inv(double(Psi_sol_global(s_i)))-eye(3)/chi_m;
    N_curr_local = double(rot_Matrix(s_i)).'*N_curr_global*double(rot_Matrix(s_i));
    N_11_vec(i) = double(N_curr_local(1,1));
    N_22_vec(i) = double(N_curr_local(2,2));
    N_33_vec(i) = double(N_curr_local(3,3));
    
end
%%
plot(s_plot,N_11_vec,s_plot,N_22_vec,s_plot,N_33_vec)
legend('N_{11}(s)','N_{22}(s)','N_{33}(s)')
grid on
axis tight
xlabel('s')



%% Plotting of M(s) in Local Coordinates
% close all
% figure(1)
% ezplot(inv(Psi_sol_local_plot(1,1))-eye(3)/chi_m,[s_min,s_max])
% hold on
% ezplot(inv(Psi_sol_local_plot(2,2))-eye(3)/chi_m,[s_min,s_max])
% ezplot(inv(Psi_sol_local_plot(3,3))-eye(3)/chi_m,[s_min,s_max])
% hold off
% legend('N_{1,1}(s)','N_{2,2}(s)','N_{3,3}(s)')
% title({'DemagFactor N(s) in local coordinates',...
%     ['Helix: H' num2str(n_helix)]})
% axis tight
% grid on
%
%
