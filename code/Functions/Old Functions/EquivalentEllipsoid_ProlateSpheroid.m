function s = EquivalentEllipsoid_ProlateSpheroid(N_11)

% s = b/a and b = c
% There are no constraints between a and b = c, the demag factor used here
% will be the one corresponding to the unequal axis a.
% N_11 corresponds to the demagnetizing factor in direction of 'a'
% N_22 = N_33 corresponds to the demagnetizing factor in direction ob 'b'
% N_33 = 1 - N_22 - N_11 corresponds to the demagnetizing factor of 'c'

%Source: Demagnetizing Factors of the General Ellipsoid, J.A. Osborn, Noval
%Ordenance Laboratory, Washington, D.C., March 24, 1945


% m = sym('m');
% assume(m, 'real');
% assume(m>= 0);

%L = @(m) 1/(m^2-1)*(m/(2*sqrt(m^2-1))*log((m+sqrt(m^2-1))/(m-sqrt(m^2-1)))-1);
%s = 1/vpasolve(L(m) == N_11, m);

%m_0 = 1;
%s = 1/fzero(L-N_11,m_0);

% L(m) is a function with D(L) = [0,inf) and W(L) = (0,1] and is
% monotonically decaying.

problem.objective = @(m) 1/(m^2-1)*(m/(2*sqrt(m^2-1))*log((m+sqrt(m^2-1))/(m-sqrt(m^2-1)))-1)- N_11;
problem.x0 = 2;
problem.solver = 'fzero'; % a required part of the structure
problem.options = optimset(@fzero); % default options
s = 1/fzero(problem);

end
