function  N = DemagCylinder(Beta,dir)
%   DemagCylindar calculates the demagfactors for a cylindar
%   Beta specifies the Diamter/Length ratio
%   dir = 'Axial' or 'Radial' specifes the factor desired

N = zeros(length(Beta));

for i = 1:length(Beta)
    
    Eta =1/Beta(i);
    Keta = 1/sqrt(1+Eta^2);
    [K,E] = ellipke(Keta^2);
    if strcmp(dir, 'Axial')
        N(i) = 1+4/(3*pi*Eta)*(1-1/Keta*((1-Eta^2)*E+Eta^2*K));
    else
        N(i) = -2/(3*pi*Eta)*(1-1/Keta*((1-Eta^2)*E+Eta^2*K));
    end
    
end
end
