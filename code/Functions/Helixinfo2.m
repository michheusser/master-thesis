function [h_helix, D_helix, k_helix, d_helix] = Helixinfo2(n_helix)

d_nom = 3e-6; %m
D_nom = 10e-6; %m
h_nom = 10e-6; %m

k_helix = 3; %1
d_helix = d_nom; %m
h_min = d_helix; %m
h_max = 32.97e-6; %m

v_nom = (pi*d_nom^2/4)*k_helix*sqrt(h_nom^2 + (pi*D_nom)^2); %m^3

steps = 10;

h_helix = round((h_min + (n_helix-1)*(h_max-h_min)/(steps-1))*1e8)*1e-8;
D_helix = round(((1/pi)*sqrt((round((4*v_nom/(k_helix*pi*d_nom^2))*1e8)*1e-8)^2-h_helix.^2))*1e8)*1e-8;



