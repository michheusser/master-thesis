close all
clear all
clc


%%

d_helix = 3e-6;
k_helix = 3;


DF = load('DemagFactors.mat')
names = fieldnames(DF.DemagFactors)
Helix_vec = [1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 4 5 6 7 8 9 10]';

for i = 16 : 16
    N = DF.DemagFactors.(['S' num2str(i)]).N_simulationM;
    
    [V,D] = eig(N);
    N_x = D(1,1);
    N_y = D(2,2);
    
    
    
    n_helix = Helix_vec(i);
    [h_helix, D_helix] = Helixinfo(n_helix);
    vol = (pi*d_helix^2/4)*k_helix*sqrt(h_helix^2 + (pi*D_helix)^2);
    [a, b, c] = EquivalentEllipsoid(N_x,N_y,vol)
[x, y, z] = ellipsoid(0,0,0,a,b,c,30);
x_rot = x;
y_rot = y;
z_rot = z;

[m,n] = size(x);

for ii = 1 : m
   
    for jj = 1: n
       
       v_old = [x(ii,jj);y(ii,jj);z(ii,jj)];
       v_new = V'*v_old;
       x_rot(ii,jj) = v_new(1);
       y_rot(ii,jj) = v_new(2);  
       z_rot(ii,jj) = v_new(3);
        
       
    end
    
    
end

end

%%

x_range = [-1e-5 1e-5];
y_range = x_range;
z_range = x_range;

close all
subplot(1,2,1)
S = surfl(x, y, z);
xlim(x_range)
ylim(y_range)
zlim(z_range)

subplot(1,2,2)
S = surfl(x_rot, y_rot, z_rot);
xlim(x_range)
ylim(y_range)
zlim(z_range)



