

syms cphi sphi ctheta stheta calpha salpha
sphi = sqrt(1-cphi^2);
stheta = sqrt(1-ctheta^2);
salpha = sphi/stheta;
calpha = sqrt(1-salpha^2);
theta = acos(ctheta);

N_11 = 0.333333
N_22 = 0.333333

P = vpasolve([cphi*ctheta/(stheta^3*salpha^2)*(ellipticF(theta,salpha)-ellipticE(theta,salpha))== N_11, cphi*ctheta/(sphi^3*calpha^2)*(stheta*cphi/ctheta - ellipticE(theta,salpha)) == N_22], [cphi, ctheta])

P.cphi
P.ctheta