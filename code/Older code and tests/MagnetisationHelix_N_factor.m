close all
clear all
clc

%% Function Inputs

N = 10;
M_sol_global_handle_cell = cell(3,1);
M_sol_local_handle_cell = cell(3,1);
M_sol_av_cell = cell(3,1);

chi_m = 1000;

integration_method = 'auto';
shape_function = 'fourier';
n_helix = 3;
B_app_matrix = eye(3);

for i = 1 : 3
    
    B_app = B_app_matrix(:,i);
    [M_sol_global_handle, M_sol_local_handle, M_sol_av] = Magnetisation_Helix(n_helix,B_app,chi_m,N,shape_function, integration_method);
    M_sol_global_handle_cell{i,1} = M_sol_global_handle;
    M_sol_local_handle_cell{i,1} = M_sol_local_handle;
    M_sol_av_cell{i,1} = M_sol_av;
    
end

save('MagnetizationH3Bxyz.mat')

%% Plot all

close all
f1 = figure('Name',['Magnetization M(s) for N = ' num2str(N)]);
a1 = axes('Parent',f1);
h = zeros(3,3);
style_cell = {':';'-.';'--'};
color_cell = {'b';'g';'r'};

[~,~,k_helix,~] = Helixinfo(n_helix);
s_min = 0;
s_max = 2*pi*k_helix;
legend_cell = cell(3,1);

hold(a1,'on')
for i = 1 : 3
    M_sol_local_sym = sym(M_sol_local_handle_cell{i,1});
    
    for j = 1 : 3
        M_sol_curr = matlabFunction(M_sol_local_sym(j));
        h(i,j) = ezplot(a1,M_sol_curr,[s_min,s_max]);
        set(h(i,j),'LineStyle',style_cell{i,1})
        set(h(i,j),'Color',color_cell{j,1})
        legend_cell{(i-1)*3+j,1} = ['M_{' num2str(j) '}(s)  - (B in direction ' num2str(i) ' ) N = ' num2str(N)];
    end
    
end
hold(a1,'off')
legend(a1,legend_cell)
title(a1,{'','Magnetization M(s) in local coordinates',...
    ['Method: ' shape_function],...
    ['Helix: H' num2str(n_helix)]})
axis tight
grid on

%%

f2 = figure('Name',['Magnetization M(s) for N = ' num2str(N)]);


h = zeros(3,3);
color_cell = {'b';'g';'r'};

[~,~,k_helix,~] = Helixinfo(n_helix);
s_min = 0;
s_max = 2*pi*k_helix;
legend_cell = cell(3,1);


for i = 1 : 3
    a(i) = subplot(3,1,i,'Parent',f2);
    hold(a(i),'on')
    
    for j = 1 : 3
        M_sol_local_sym = sym(M_sol_local_handle_cell{j,1});
        M_sol_curr = matlabFunction(M_sol_local_sym(i));
        h(i,j) = ezplot(a(i),M_sol_curr,[s_min,s_max]);
        set(h(i,j),'Color',color_cell{j,1})
 %       legend_cell{(i-1)*3+j,1} = ['M_{' num2str(j) '}(s)  - (B in direction ' num2str(i) ' ) N = ' num2str(N)];
    end
    axis tight
grid on
hold(a(i),'off')    
end

%legend(a1,legend_cell)
%title(a1,{'','Magnetization M(s) in local coordinates',...
    %['Method: ' shape_function],...
    %['Helix: H' num2str(n_helix)]})



%% Calculation of N-Factor

mu_0 = 1.25663706e-6;
H_app_matrix = B_app_matrix./mu_0;

N_s = 20;
s_i_vec = s_min + (1:N_s)/(N_s+1)*(s_max-s_min);
N_matrix_global_line_cell = cell(length(s_i_vec),1);
N_matrix_local_line_cell = cell(length(s_i_vec),1);



for i = 1 : length(s_i_vec)
    
    s_i = s_i_vec(i);
    M_x_global_handle = M_sol_global_handle_cell{1,1}; 
    M_y_global_handle = M_sol_global_handle_cell{2,1};
    M_z_global_handle = M_sol_global_handle_cell{3,1};
    M_global_matrix = [M_x_global_handle(s_i) M_y_global_handle(s_i) M_z_global_handle(s_i)];
    N_matrix_global_curr = H_app_matrix/M_global_matrix - eye(3)/chi_m;

end




