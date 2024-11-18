clear; close all;
dbstop if error;

addpath(".\src\")

%% Parameters for composite time integration
scheme = 'M-PF';        % select "M-PF" or "PadePF"
nSubStep = 4;           % max. 6 for "M-PF"
rhoInfty = 0;           % value between 0 and 1

%% Parameters of 3-DOF system
tmax = 100;        % simulation time (seconds)
dt = 0.14;          % time step size
% Spring constant
k1 = 10^7;
k2 = 1;
% Mass
m1 = 0;
m2 = 1;
m3 = 1;

np = 2;
M = [m2 0; 0, m3];          % Mass matrix
C = zeros(2);               % Damping matrix
K = [k1+k2, -k2; -k2, k2];  % Stiffness matrix
F0 = [1; 0]*k1;             % Force vector
u0 = [0; 0];                % Initial displacement
v0 = [0; 0];                % Initial velocity
BC_Accl = [];

% Save results for particular degrees of freedom
pDOF = [1;2];

% Excitation function
amp = 1;                            % amplitude
omega_ex = 1.2;                     % excitation frequency
fHist = @(t) amp*(sin(omega_ex*(t)));      % excitation signal


%% Reference solution
tp = (0:0.02:tmax);
[uRef, vRef, aRef] = ThreeDOFsRefSln(M,K,F0,tp);
R1Ref = k1*(fHist(tp)-uRef(1,:));

%% Time stepping solution
ns = floor(tmax/dt) + 1;    % number of time steps

% store user-parameters
user_params = [nSubStep, rhoInfty];

switch scheme
    
    % scheme using a single multiple root
    case 'M-PF'
        [dsp, vel, acc] = TimeSolverTRhoPF(user_params,fHist, ...
                                        ns,dt,K,M,C,F0,u0,v0,pDOF);

    % scheme using distinct roots
    case 'PadePF'
        [dsp, vel, acc] = TimeSolverPF(user_params,fHist,...
                                        ns,dt,K,M,C,F0,u0,v0,pDOF);

end
tn = (0:ns-1)*dt;
R1 = k1*((fHist(tn))'-dsp(:,1));


%% Plot results

figure(1)
plot(tp,uRef(1,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,dsp(:,1), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Displacement')
title('Displacement of Mass 2');

figure(2)
plot(tp,uRef(2,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,dsp(:,2), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Displacement')
title('Displacement of Mass 3');

figure(3)
plot(tp,vRef(1,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,vel(:,1), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Velocity')
title('Velocity of Mass 2');

figure(4)
plot(tp,vRef(2,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,vel(:,2), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Velocity')
title('Velocity of Mass 3');

figure(5)
plot(tp,aRef(1,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,acc(:,1), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Acceleration')
title('Acceleration of Mass 2');

figure(6)
plot(tp,aRef(2,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,acc(:,2), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Acceleration')
title('Acceleration of Mass 3');

figure(7)
plot(tp,R1Ref, '-r',"DisplayName",'Reference')
hold on
plot(tn,R1, '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Force')
title('Reaction Force');


%% Reference solution
function [uRef, vRef, aRef] = ThreeDOFsRefSln(M,K,F0,tp)
[Vec,D] = eig(full(K));
Mg = Vec'*M*Vec;
Kg = Vec'*K*Vec;
Fg = Vec'*F0;
o1 = sqrt(D(1,1));
o2 = sqrt(D(2,2));
Uex = @(o,t) 1/(o*o-1.2^2)*(sin(1.2*t) - 1.2/o*sin(o*t));
Uref = @(o,t) 1/(o*o-1.2^2)*(sin(1.2*t)       );
Vex = @(o,t) 1.2/(o*o-1.2^2)*(cos(1.2*t) - cos(o*t));
Vref = @(o,t) 1.2/(o*o-1.2^2)*(cos(1.2*t)     );
Aex = @(o,t) 1.2/(o*o-1.2^2)*(-1.2*sin(1.2*t) + o*sin(o*t));
Aref = @(o,t) 1.2/(o*o-1.2^2)*(-1.2*sin(1.2*t)     );
U  = [Fg(1)*Uex(o1,tp); Fg(2)*Uref(o2,tp)];
uRef = Vec*U;
V  = [Fg(1)*Vex(o1,tp); Fg(2)*Vref(o2,tp)];
vRef = Vec*V;
A  = [Fg(1)*Aex(o1,tp); Fg(2)*Aref(o2,tp)];
aRef = Vec*A;
end

