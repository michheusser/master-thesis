close all
clear all
clc

%% Raw

Helix_Data_Raw_2 = load('NiHelix_k2/Helix_Data_Raw');
Helix_Data_Raw_Struct_2 = Helix_Data_Raw_2.Helix_Data_Raw;
M_ax_ax_2 = Helix_Data_Raw_Struct_2.M_ax_ax;
M_rad_rad_2 = Helix_Data_Raw_Struct_2.M_rad_rad;

Helix_Data_Raw_3 = load('NiHelix_k3/Helix_Data_Raw');
Helix_Data_Raw_Struct_3 = Helix_Data_Raw_3.Helix_Data_Raw;
M_ax_ax_3 = Helix_Data_Raw_Struct_3.M_ax_ax;
M_rad_rad_3 = Helix_Data_Raw_Struct_3.M_rad_rad;

Helix_Data_Raw_4 = load('NiHelix_k4/Helix_Data_Raw');
Helix_Data_Raw_Struct_4 = Helix_Data_Raw_4.Helix_Data_Raw;
M_ax_ax_4 = Helix_Data_Raw_Struct_4.M_ax_ax;
M_rad_rad_4 = Helix_Data_Raw_Struct_4.M_rad_rad;

Helix_Data_Raw_7 = load('NiHelix_k7/Helix_Data_Raw');
Helix_Data_Raw_Struct_7 = Helix_Data_Raw_7.Helix_Data_Raw;
B_app = Helix_Data_Raw_Struct_7.B_app_ax;
M_ax_ax_7 = Helix_Data_Raw_Struct_7.M_ax_ax;
M_rad_rad_7 = Helix_Data_Raw_Struct_7.M_rad_rad;

plot(B_app,M_ax_ax_2,'b-',B_app,M_rad_rad_2,'b--',...
    B_app,M_ax_ax_3,'g-',B_app,M_rad_rad_3,'g--',...
    B_app,M_ax_ax_4,'r-',B_app,M_rad_rad_4,'r--',...
    B_app,M_ax_ax_7,'k-',B_app,M_rad_rad_7,'k--')

grid on
axis tight
xlabel('Applied Field H [Oe]')
ylabel('Magnetic moment m [memu]')
legend('k = 2, Axial','k = 2, Radial',...
    'k = 3, Axial','k = 3, Radial',...
    'k = 4, Axial','k = 4, Radial',...
    'k = 7, Axial','k = 7, Radial')
title('VSM Measurements for various helices (Magnetic Moment)')


%% Not Normalized

Helix_Data_2 = load('NiHelix_k2/Helix_Data');
Helix_Data_Struct_2 = Helix_Data_2.Helix_Data;
M_ax_ax_2 = Helix_Data_Struct_2.M_ax_ax;
M_rad_rad_2 = Helix_Data_Struct_2.M_rad_rad;

Helix_Data_3 = load('NiHelix_k3/Helix_Data');
Helix_Data_Struct_3 = Helix_Data_3.Helix_Data;
M_ax_ax_3 = Helix_Data_Struct_3.M_ax_ax;
M_rad_rad_3 = Helix_Data_Struct_3.M_rad_rad;

Helix_Data_4 = load('NiHelix_k4/Helix_Data');
Helix_Data_Struct_4 = Helix_Data_4.Helix_Data;
M_ax_ax_4 = Helix_Data_Struct_4.M_ax_ax;
M_rad_rad_4 = Helix_Data_Struct_4.M_rad_rad;

Helix_Data_7 = load('NiHelix_k7/Helix_Data');
Helix_Data_Struct_7 = Helix_Data_7.Helix_Data;
B_app = Helix_Data_Struct_7.B_app_ax;
M_ax_ax_7 = Helix_Data_Struct_7.M_ax_ax;
M_rad_rad_7 = Helix_Data_Struct_7.M_rad_rad;

plot(B_app,M_ax_ax_2,'b-',B_app,M_rad_rad_2,'b--',...
    B_app,M_ax_ax_3,'g-',B_app,M_rad_rad_3,'g--',...
    B_app,M_ax_ax_4,'r-',B_app,M_rad_rad_4,'r--',...
    B_app,M_ax_ax_7,'k-',B_app,M_rad_rad_7,'k--')

grid on
axis tight
xlabel('Applied Field B [T]')
ylabel('Magnetization M [A/m]')
legend('k = 2, Axial','k = 2, Radial',...
    'k = 3, Axial','k = 3, Radial',...
    'k = 4, Axial','k = 4, Radial',...
    'k = 7, Axial','k = 7, Radial')
title('VSM Measurements for various helices (Magnetization)')


%% Normalized

Helix_Data_2 = load('NiHelix_k2/Helix_Data_Normalized');
Helix_Data_Struct_2 = Helix_Data_2.Helix_Data_Normalized;
M_ax_ax_2 = Helix_Data_Struct_2.M_ax_ax;
M_rad_rad_2 = Helix_Data_Struct_2.M_rad_rad;

Helix_Data_3 = load('NiHelix_k3/Helix_Data_Normalized');
Helix_Data_Struct_3 = Helix_Data_3.Helix_Data_Normalized;
M_ax_ax_3 = Helix_Data_Struct_3.M_ax_ax;
M_rad_rad_3 = Helix_Data_Struct_3.M_rad_rad;

Helix_Data_4 = load('NiHelix_k4/Helix_Data_Normalized');
Helix_Data_Struct_4 = Helix_Data_4.Helix_Data_Normalized;
M_ax_ax_4 = Helix_Data_Struct_4.M_ax_ax;
M_rad_rad_4 = Helix_Data_Struct_4.M_rad_rad;

Helix_Data_7 = load('NiHelix_k7/Helix_Data_Normalized');
Helix_Data_Struct_7 = Helix_Data_7.Helix_Data_Normalized;
B_app = Helix_Data_Struct_7.B_app_ax;
M_ax_ax_7 = Helix_Data_Struct_7.M_ax_ax;
M_rad_rad_7 = Helix_Data_Struct_7.M_rad_rad;

plot(B_app,M_ax_ax_2,'b-',B_app,M_rad_rad_2,'b--',...
    B_app,M_ax_ax_3,'g-',B_app,M_rad_rad_3,'g--',...
    B_app,M_ax_ax_4,'r-',B_app,M_rad_rad_4,'r--',...
    B_app,M_ax_ax_7,'k-',B_app,M_rad_rad_7,'k--')

grid on
axis tight
xlabel('Applied Field B [T]')
ylabel('Normalized Magnetization M/M_s [1]')
legend('k = 2, Axial','k = 2, Radial',...
    'k = 3, Axial','k = 3, Radial',...
    'k = 4, Axial','k = 4, Radial',...
    'k = 7, Axial','k = 7, Radial')
title('VSM Measurements for various helices (Normalized Magnetization)')

