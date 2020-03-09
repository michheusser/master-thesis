function [m, n] = EquivalentEllipsoid(N_11,N_22)

% m and n are the relations of the ellipsoid axes. Where m = b/a and n = c/a
% and a >= b >= c.
% N_11, N_22 are the demagnetizing factors of any of the two easy axes of
% the demagnetizing tensor. The function calculates the third demag factor
% N_33 = 1 - N_22 - N_11 and calculates m and n under the assumption that
% a >=b >= c, where L corresponds to a (the biggest), M corresponds to b
% (medium size) and N corresponds to c (the smallest)

%Source Demagnetizing Factors of the General Ellipsoid, J.A. Osborn, Noval
%Ordenance Laboratory, Washington, D.C., March 24, 1945


N_33 = 1 - N_11 - N_22;

N = [N_11; N_22; N_33]; %N_11 biggest, N_33 smallest
N = sort(N,'descend');

%%
% syms cphi sphi ctheta stheta salpha calpha theta k 
% sphi = sqrt(1-cphi^2);
% stheta = sqrt(1-ctheta^2);
% salpha = sphi/stheta;
% calpha = sqrt(1-salpha^2);
% theta = acos(ctheta);
% k = salpha;
% 
% L(cphi,ctheta) = cphi*ctheta/(stheta^3*salpha^2)*(ellipticF(theta,k^2)-ellipticE(theta,k^2));
% M(cphi,ctheta) = cphi*ctheta/(stheta^3*salpha^2*calpha^2)*(ellipticE(theta,k^2)-calpha^2*ellipticF(theta,k^2)-salpha^2*stheta*ctheta/cphi);
% %N(cphi,ctheta) = cphi*ctheta/(sphi^3*calpha^2)*(stheta*cphi/ctheta - ellipticE(theta,k^2));
% 
% P = vpasolve([L(cphi,ctheta) == N_11, M(cphi,ctheta) == N_22],[cphi,ctheta]);
% %P = vpasolve([L(cphi,ctheta) == N_11, N(cphi,ctheta) == N_33],[cphi,ctheta]);
% %P = vpasolve([M(cphi,ctheta) == N_22, N(cphi,ctheta) == N_33],[cphi,ctheta]);
% 
% m = P.cphi; % m = b/a
% n = P.ctheta; % n = c/a


%L(m(1),m(2)) = m(1)*m(2)/(sqrt(1-m(2)^2)^3*(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)*(ellipticF(acos(m(2)),(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)-ellipticE(acos(m(2)),(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2));
%M(m(1),m(2)) = m(1)*m(2)/(sqrt(1-m(2)^2)^3*(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2*sqrt(1-(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)^2)*(ellipticE(acos(m(2)),(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)-sqrt(1-(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)^2*ellipticF(acos(m(2)),(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)-(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2*sqrt(1-m(2)^2)*m(2)/m(1));


%m_0 = EquivalentEllipsoid_ProlateSpheroid(N_11)

problem.objective = @(m) [m(1)*m(2)/(sqrt(1-m(2)^2)^3*(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)*(ellipticF(acos(m(2)),(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)-ellipticE(acos(m(2)),(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)) - N_11; m(1)*m(2)/(sqrt(1-m(2)^2)^3*(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2*sqrt(1-(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)^2)*(ellipticE(acos(m(2)),(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)-sqrt(1-(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)^2*ellipticF(acos(m(2)),(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2)-(sqrt(1-m(1)^2)/sqrt(1-m(2)^2))^2*sqrt(1-m(2)^2)*m(2)/m(1)) - N_22] 
problem.x0 = [0.1;0.1];
problem.solver = 'fsolve'; % a required part of the structure
problem.options = optimset(@fsolve); % default options
result = fsolve(problem);
m = result(1);
n = result(2);



end

