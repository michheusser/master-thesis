function N = DMG_Rectangular_Prism(a,b,c)
%WARNING: a,b,c are the half-length of the sides of the rectangle in x,y
% and z axes, respectively
% N is the demagnetization factor in z direction for an applied b-field in
% z direction. 

%Aharoni, 1998

d = sqrt(a^2+b^2+c^2);
alpha = sqrt(a^2+b^2);
beta = sqrt(b^2+c^2);
gamma = sqrt(c^2+a^2);


N = ((b^2-c^2)/(2*b*c)*log((d-a)/(d+a)) + ...
    (a^2-c^2)/(2*a*c)*log((d-b)/(d+b)) + ...
    b/(2*c)*log((alpha+a)/(alpha-a)) + ...
    a/(2*c)*log((alpha+b)/(alpha-b)) + ...
    c/(2*a)*log((beta-b)/(beta+b)) + ...
    c/(2*b)*log((gamma-a)/(gamma+a)) + ...
    2*atan((a*b)/(c*d)) + ...
    (a^3 + b^3 - 2*c^3)/(3*a*b*c) + ...
    (a^2 + b^2 - 2*c^2)/(3*a*b*c)*d + ...
    c/(a*b)*(gamma+beta) - ...
    (alpha^3 + beta^3 + gamma^3)/(3*a*b*c))/pi;



end
