function [N_x,N_y,N_z] = Demagfactor_Ellipsoid_General(tau_a,tau_b)

% tau_a := c/a, tau_b := c/b
% a -> N_x, b -> N_y, c -> N_z
if(tau_a == 1 && tau_b ~= 1) % (a = c)
    
    N_y = Demagfactor_Ellipsoid_Unidirectional(tau_a/tau_b,1/tau_b);
    N_z = 0.5*(1-N_y);
    N_x = N_z;
    
elseif(tau_a ~= 1 && tau_b == 1) %(b = c)
    
    N_x = Demagfactor_Ellipsoid_Unidirectional(1/tau_a,tau_b/tau_a);
    N_y = 0.5*(1-N_x);
    N_z = N_y;
    
elseif(tau_a == tau_b) %(a = b)
    
    N_z = Demagfactor_Ellipsoid_Unidirectional(tau_a,tau_b);
    N_x = 0.5*(1-N_z);
    N_y = N_x;
    
else
    
    N_z = Demagfactor_Ellipsoid_Unidirectional(tau_a,tau_b);
    N_x = Demagfactor_Ellipsoid_Unidirectional(1/tau_a,tau_b/tau_a);
    N_y = Demagfactor_Ellipsoid_Unidirectional(tau_a/tau_b,1/tau_b);
    
end
