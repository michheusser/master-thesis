function N = Demagfactor_Ellipsoid_Unidirectional(t_a,t_b)

% tau_a := c/a, tau_b := c/b
% a -> N_x, b -> N_y, c -> N_z
% Source: The equivalent ellipsoid of a magnetized body, M Beleggia et al.,
%  2006

if(t_a == t_b) % For the cases where t_a = 1 (c = a), t_b = 1 (c = b) or t_a == t_b (a = b)
    t_e = t_a;
    
    if(t_e == 1)
        N = 1/3;
    else
        N = 1/(1-t_e^2)*(1-(t_e*acos(t_e))/(sqrt(1-t_e^2)));
    end
    
else % For all other cases
    k = asin(sqrt(1-t_a^(-2)));
    m = (1 - t_b^(-2))/(1 - t_a^(-2));
    N = 1/(t_a*t_b)*(ellipticF(k,m)-ellipticE(k,m))/(m*sin(k)^3);
end

