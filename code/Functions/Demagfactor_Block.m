function N = Demagfactor_Block(a,b,c)
%a (x-axis), b (y-axis) c (z-axis)
N_33_analytical = DMG_Rectangular_Prism(a,b,c);
N_11_analytical = DMG_Rectangular_Prism(c,a,b);
N_22_analytical = DMG_Rectangular_Prism(b,c,a);

N = diag([N_11_analytical, N_22_analytical, N_33_analytical]);

end

