close all
clear all
clc

%% Function Inputs

N_vec = [3 5 10 15];
M_sol_global_handle_cell = cell(length(N_vec),1);
M_sol_local_handle_cell = cell(length(N_vec),1);
M_sol_av_cell = cell(length(N_vec),1);

for i = 1 : length(N_vec)
    N = N_vec(i);
    B_app = [0 0 1]';
    n_helix = 3;
    chi_m = 1000;
    shape_function = 'polynomials';
    [M_sol_global_handle, M_sol_local_handle, M_sol_av] = Magnetisation_Helix(n_helix,B_app,chi_m,N,shape_function);
    M_sol_global_handle_cell{i,1} = M_sol_global_handle;
    M_sol_local_handle_cell{i,1} = M_sol_local_handle;
    M_sol_av_cell{i,1} = M_sol_av;
end

%%
close all
f1 = figure('Name','Comparison between approximations of Magnetization M(s)');
a1 = axes('Parent',f1);
h = zeros(length(N_vec),3);
style_cell = {':';'-.';'--';'-'};
color_cell = {'b';'g';'r'};

[~,~,k_helix,~] = Helixinfo(n_helix);
s_min = 0;
s_max = 2*pi*k_helix;
legend_cell = cell(length(N_vec)*3,1);

hold(a1,'on')
for i = 1 : length(N_vec)
    M_sol_local_sym = sym(M_sol_local_handle_cell{i,1});
    
    for j = 1 : 3
        M_sol_curr = matlabFunction(M_sol_local_sym(j));
     h(i,j) = ezplot(a1,M_sol_curr,[s_min,s_max]);
    set(h(i,j),'LineStyle',style_cell{i,1})
    set(h(i,j),'Color',color_cell{j,1})
    legend_cell{(i-1)*3+j,1} = ['M_{' num2str(j) '}(s)  - N = ' num2str(N_vec(i))];
    end
    
end
hold(a1,'off')
legend(a1,legend_cell)
title(a1,{'','Magnetization M(s) in local coordinates',...
     ['Applied field: B = ( ' num2str(B_app(1)) ', ' num2str(B_app(2)) ', ' num2str(B_app(3)) ')^T [T]'],...
     ['Method: ' shape_function],...
     ['Helix: H' num2str(n_helix)]})
 axis tight
 grid on
 
 %%