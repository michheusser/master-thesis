function N = Demagfactor_Ellipsoid_General_2(m,n)



N_11 = m/(2*(m^2-1))*(m-(1/(2*sqrt(m^2-1)))*log((m+sqrt(m^2-1))/(m-sqrt(m^2-1))));
N_22 = N_11;
N_33 = 1/(m^2-1)*(m/(2*sqrt(m^2-1))*log((m+sqrt(m^2-1))/(m-sqrt(m^2-1)))-1);

N = diag([N_11, N_22 ,N_33]);