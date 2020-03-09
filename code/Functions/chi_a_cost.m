function chi_a_cost_val = chi_a_cost( chi_m,  N_simulation , chi_a)

chi_a_cost_val = zeros(length(chi_m),1);

for i = 1 : length(chi_m)
chi_a_matrix = eye(3)*((eye(3)/chi_m(i)) + N_simulation)^(-1) - chi_a;
chi_a_cost_val(i) = (chi_a_matrix(1,1)^2 + chi_a_matrix(2,2)^2 + 0*chi_a_matrix(3,2)^2 + 0*chi_a_matrix(2,3)^2);

end

end

