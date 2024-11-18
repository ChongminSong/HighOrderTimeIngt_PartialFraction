clear; close all;
dbstop if error;

addpath(".\src\")

%% Parameters for composite time integration
scheme = 'PadePF';        % select "M-PF" or "PadePF"
nSubStep = 4;           % max. 6 for "M-PF"
rhoInfty = 1;           % value between 0 and 1

%% Parameters of SDOF system
tmax = 10;
    
omega = 2*pi;               % natural (angular) frequency [rad/s]
zeta = 0;                   % damping ratio [%]
m = 1;                      % mass [kg]

% Derived parameters
k = omega^2*m;      % stiffness [N/m]
c = 2*zeta*omega*m;	% damping coefficient [Ns/m]
f = omega/(2*pi);   % natural frequency [Hz]
dt = f/40;           % time increment [s]

% Initial conditions
u0 = 2;                     % initial displacement
v0 = pi/3;                  % initial velocity

% Loading
omega_ex1 = 2*sqrt(5)/5;	% excitation frequency (cos-term) [rad/s]
omega_ex2 = 2*sqrt(10);     % excitation frequency (sin-term) [rad/s]
a1 = 10;                    % force amplitude (cos-term) [N]
a2 = 70;                    % force amplitude (sin-term) [N]
fHist = @(t) a1*(cos(omega_ex1*(t)))+a2*(sin(omega_ex2*(t)));      
F0 = 1;
BC_Accl = [];

% Save results for particular degrees of freedom
pDOF = 1;     	% degrees of freedom whose results are saved
nDOF = 1;           % number of degrees of freedom

%% Reference solution
C1 = u0 - (a1/k)/(1-(omega_ex1/omega)^2);
C2 = v0/omega - a2/k*(omega_ex2/omega)/(1-(omega_ex2/omega)^2);
C3 = (a1/k)/(1-(omega_ex1/omega)^2);
C4 = (a2/k)/(1-(omega_ex2/omega)^2);
u_exact = @(t) C1*cos(omega*t)     + ...
               C2*sin(omega*t)     + ...
               C3*cos(omega_ex1*t) + ...
               C4*sin(omega_ex2*t);
v_exact = @(t) -omega*C1*sin(omega*t)         + ...
                omega*C2*cos(omega*t)         + ...
               -omega_ex1*C3*sin(omega_ex1*t) + ...
                omega_ex2*C4*cos(omega_ex2*t);
a_exact = @(t) -omega^2*C1*cos(omega*t)         + ...
               -omega^2*C2*sin(omega*t)         + ...
               -omega_ex1^2*C3*cos(omega_ex1*t) + ...
               -omega_ex2^2*C4*sin(omega_ex2*t);

tp = (0:0.02:tmax);
uRef = u_exact(tp);
vRef = v_exact(tp);
aRef = a_exact(tp);

%% Time stepping solution
ns = floor(tmax/dt) + 1;    % number of time steps

% store user-parameters
user_params = [nSubStep, rhoInfty];

switch scheme
    
    % scheme using a single multiple root
    case 'M-PF'
        [dsp, vel, acc] = TimeSolverTRhoPF(user_params,fHist, ...
                                        ns,dt,k,m,c,F0,u0,v0,pDOF);

    % scheme using distinct roots
    case 'PadePF'
        [dsp, vel, acc] = TimeSolverPF(user_params,fHist,...
                                        ns,dt,k,m,c,F0,u0,v0,pDOF);

end
tn = (0:ns-1)*dt;

%% Plot results

figure(1)
plot(tp,uRef(1,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,dsp(:,1), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Displacement')
title('Displacement');

figure(3)
plot(tp,vRef(1,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,vel(:,1), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Velocity')
title('Velocity');

figure(5)
plot(tp,aRef(1,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,acc(:,1), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Acceleration')
title('Acceleration');

