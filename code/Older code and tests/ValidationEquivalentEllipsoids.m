close all
clear all
clc
%%

%m = b/a
%n = c/a

%m_vec = [0: 0.1 : 1]';
m_vec = 0.5;
n_vec = m_vec;

L_vec = zeros(length(m_vec),1);
M_vec = zeros(length(m_vec),1);

for i = 1 : length(m_vec)
    
    m = m_vec(i);
    n = n_vec(i);
    
    
    cphi = m;
    ctheta = n;
    sphi = sqrt(1-cphi^2);
    stheta = sqrt(1-ctheta^2);
    salpha = sphi/stheta
    calpha = sqrt(1-salpha^2);
    theta = acos(ctheta);
    k = salpha;
    
    L = cphi*ctheta/(stheta^3*salpha^2)*(ellipticF(theta,k^2)-ellipticE(theta,k^2));
    M = cphi*ctheta/(stheta^3*salpha^2*calpha^2)*(ellipticE(theta,k^2)-calpha^2*ellipticF(theta,k^2)-salpha^2*stheta*ctheta/cphi);
    
    %L = m*n/(sqrt(1-n^2)^3*(sqrt(1-m^2)/sqrt(1-n^2))^2)*(ellipticF(acos(n),(sqrt(1-m^2)/sqrt(1-n^2))^2)-ellipticE(acos(n),(sqrt(1-m^2)/sqrt(1-n^2))^2));
    %M = m*n/(sqrt(1-n^2)^3*(sqrt(1-m^2)/sqrt(1-n^2))^2*sqrt(1-(sqrt(1-m^2)/sqrt(1-n^2))^2)^2)*(ellipticE(acos(n),(sqrt(1-m^2)/sqrt(1-n^2))^2)-sqrt(1-(sqrt(1-m^2)/sqrt(1-n^2))^2)^2*ellipticF(acos(n),(sqrt(1-m^2)/sqrt(1-n^2))^2)-(sqrt(1-m^2)/sqrt(1-n^2))^2*sqrt(1-n^2)*n/m)
    
    
    L_vec(i) = L;
    M_vec(i) = M;
    
end

subplot(2,1,1)
plot(m_vec,L_vec)
subplot(2,1,2)
plot(n_vec,M_vec)

