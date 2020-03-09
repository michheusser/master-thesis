tau_a_vec = 0:0.1:2;
tau_b_vec = 0:0.1:2;

N_x_vec = zeros(length(tau_b_vec),length(tau_a_vec));
N_y_vec = zeros(length(tau_b_vec),length(tau_a_vec));
N_z_vec = zeros(length(tau_b_vec),length(tau_a_vec));

for i = 1: length(tau_a_vec)
    
    for ii = 1 : length(tau_b_vec)
        [N_x, N_y, N_z] = Demagfactor_Ellipsoid_General(tau_a_vec(i),tau_b_vec(i));
        
        N_x_vec(ii,i) = N_x;
        N_y_vec(ii,i) = N_y;
        N_z_vec(ii,i) = N_z;
        
    end
    
end

subplot(3,1,1)
surf(tau_a_vec,tau_b_vec,N_x_vec)
subplot(3,1,2)
surf(tau_a_vec,tau_b_vec,N_y_vec)
subplot(3,1,3)
surf(tau_a_vec,tau_b_vec,N_z_vec)
