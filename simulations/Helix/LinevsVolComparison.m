close all
clear all
clc

N_s = 20;

Data_Full_Line_Analytical = load('DemagMatrix_Helix_Line_Average_Analytical.mat');
Data_Full_Vol_Analytical = load('DemagMatrix_Helix_Vol_Average_Analytical.mat');

DemagFactors_Full_Line_Analytical = Data_Full_Line_Analytical.DemagMatrix_Line_Average_Analytical;
DemagFactors_Full_Vol_Analytical = Data_Full_Vol_Analytical.DemagMatrix_Vol_Average_Analytical;

n_Helix_vec = [1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 4 5 6 7 8 9 10]';

N_11_line_H1_10 = zeros(10,N_s);
N_22_line_H1_10 = zeros(10,N_s);
N_33_line_H1_10 = zeros(10,N_s);

N_11_surfav_H1_10 = zeros(10,N_s);
N_22_surfav_H1_10 = zeros(10,N_s);
N_33_surfav_H1_10 = zeros(10,N_s);

s_i_vec = [1:N_s]'/(N_s+1)*3*2*pi;

for i = 1 : length(n_Helix_vec)
    i
    n_helix_line = DemagFactors_Full_Line_Analytical.(['S' num2str(i)]).n_helix;
    n_helix_surfav = DemagFactors_Full_Vol_Analytical.(['S' num2str(i)]).n_helix;
    N_matrix_global_line_cell = DemagFactors_Full_Line_Analytical.(['S' num2str(i)]).N_matrix_global_line_cell;
    N_matrix_global_surfav_cell = DemagFactors_Full_Vol_Analytical.(['S' num2str(i)]).N_matrix_global_vol_cell;
    assert(n_helix_line==n_helix_surfav)
    
    if (mod(n_helix_line,1)==0)
        for j = 1 : N_s
            j
            s_i = s_i_vec(j);
            rot_Matrix_local = Rot_Matrix(s_i,'helix',n_helix_line);
            
            N_matrix_line = rot_Matrix_local'*N_matrix_global_line_cell{j,1}*rot_Matrix_local;
            N_11_line_H1_10(n_helix_line,j) = N_matrix_line(1,1);
            N_22_line_H1_10(n_helix_line,j) = N_matrix_line(2,2);
            N_33_line_H1_10(n_helix_line,j) = N_matrix_line(3,3);
            
            N_matrix_surfav = rot_Matrix_local'*N_matrix_global_surfav_cell{j,1}*rot_Matrix_local;
            N_11_surfav_H1_10(n_helix_line,j) = N_matrix_surfav(1,1);
            N_22_surfav_H1_10(n_helix_line,j) = N_matrix_surfav(2,2);
            N_33_surfav_H1_10(n_helix_line,j) = N_matrix_surfav(3,3);
        end
        
    end
end
%%
Helices_to_show = [1,2,3];
subplot(3,1,1)
legend_string_cell = cell(1,length(Helices_to_show)*2);
hold on
for i = 1:length(Helices_to_show)
    n_helix = Helices_to_show(i);
    hue = i/length(Helices_to_show);
    sat = 0.7;
    bright = 0.5;
    plot((1:N_s)/(N_s+1),N_11_line_H1_10(n_helix,:),'*-','Color',hsv2rgb([hue sat bright]))
    plot((1:N_s)/(N_s+1),N_11_surfav_H1_10(n_helix,:),'s--','Color',hsv2rgb([hue sat bright]))
    legend_string_cell{1,(i*2-1)} = ['H' num2str(n_helix) ': Line Point']; 
    legend_string_cell{1,(i*2)} = ['H' num2str(n_helix) ': Surface Average'];
end

hold off
grid on
xlabel('s/L (L: Arc Length)')
ylabel('N')
ylim([0 1])
title({'','Comparison between line demagnetization tensors', 'and surface average demagnetization tensors','for different helices in local coordinates','','N_{11}'})
legend(legend_string_cell,'Location','Best')

subplot(3,1,2)
hold on
for i = 1:length(Helices_to_show)
    n_helix = Helices_to_show(i);
    hue = i/length(Helices_to_show);
    sat = 0.7;
    bright = 0.5;
    plot((1:N_s)/(N_s+1),N_22_line_H1_10(n_helix,:),'*-','Color',hsv2rgb([hue sat bright]))
    plot((1:N_s)/(N_s+1),N_22_surfav_H1_10(n_helix,:),'s--','Color',hsv2rgb([hue sat bright]))
end
hold off
grid on
xlabel('s/L (L: Arc Length)')
ylabel('N')
ylim([0 1])
title({'N_{22}'})
legend(legend_string_cell,'Location','Best')

subplot(3,1,3)
hold on
for i = 1:length(Helices_to_show)
    n_helix = Helices_to_show(i);
    hue = i/length(Helices_to_show);
    sat = 0.7;
    bright = 0.5;
    plot((1:N_s)/(N_s+1),N_33_line_H1_10(n_helix,:),'*-','Color',hsv2rgb([hue sat bright]))
    plot((1:N_s)/(N_s+1),N_33_surfav_H1_10(n_helix,:),'s--','Color',hsv2rgb([hue sat bright]))
end
hold off
grid on
xlabel('s/L (L: Arc Length)')
ylabel('N')
ylim([0 1])
title({'N_{33}'})
legend(legend_string_cell,'Location','Best')