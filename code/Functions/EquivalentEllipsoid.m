function [a, b, c, theta, phi, V_std, N_std] = EquivalentEllipsoid(N,vol)

% Demagnetization Matrix: Not necessarily diagonalised
% vol: Volume of body
% a,b,c Length of semiaxes of the ellipsoid
% theta: Polar angle of ellipsoid
% phi: Azimuth angle of the ellipsoid
% V_std: Vectors of easy axes
% N_std: Diagonalised matrix

[V, N_diag] = eig(N);

N_vec = diag(N_diag);
N_std = zeros(3,3);
V_std = zeros(3,3);

[N_std(3,3), I_z] = max(N_vec);
V_std(:,3) = StdDirection(V(:,I_z));
% Ellipsoid main direction chosen

n = 1;
for i = 1 : 3
   
    if(i ~= I_z)
    N_std(n,n) = N_vec(i);
    V_std(:,n) = StdDirection(V(:,i));
    n = n+1;
    end
    
end

if(cross(V_std(:,1),V_std(:,2)).'*V_std(:,3)<0)
    
    V_std_1_tmp = V_std(:,1);
    V_std(:,1) = V_std(:,2);
    V_std(:,2) = V_std_1_tmp;

    N_std_1_tmp = N_std(1,1);
        N_std(1,1) = N_std(2,2);
    N_std(2,2) = N_std_1_tmp;
    
end

theta = 360*acos(V_std(3,3))/(2*pi);
phi = 360*acos(V_std(1,1))/(2*pi);

N_1 = N_std(1,1);
N_2 = N_std(2,2);
%phi2 = 360*acos(V_std(2,2))/(2*pi) %Should be the same as phi

% plot3([0 V_std(1,1)],[0 V_std(2,1)], [0 V_std(3,1)],...
%     [0 V_std(1,2)],[0 V_std(2,2)], [0 V_std(3,2)],...
%     [0 V_std(1,3)],[0 V_std(2,3)], [0 V_std(3,3)])
% grid on
% axis tight
% legend('V_1','V_2','V_3')
% xlabel('x')
% ylabel('y')
% zlabel('z')

tau_a = NaN;
tau_b = NaN;

if(N_1 > 0 && N_2 > 0 && N_1 + N_2 < 1)
    problem.objective = @(tau) [Demagfactor_Ellipsoid_General_Solve_x(tau)-N_1; Demagfactor_Ellipsoid_General_Solve_y(tau)-N_2];
    problem.x0 = [0.1;0.2];
    problem.solver = 'fsolve'; % a required part of the structure
    problem.options = optimoptions('fsolve','Display','off'); % default options
    result = fsolve(problem);
    tau_a = result(1);
    tau_b = result(2);
else
    
    disp('Demagnetization factors not allowed!')
    
end

c = (vol*0.75*tau_a*tau_b/pi)^(1/3);
a = c/tau_a;
b = c/tau_b;

