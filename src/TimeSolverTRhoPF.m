function [dsp,vel,acc] = TimeSolverTRhoPF(user_params,signal,ns,dt, ...
                                       K,M,C,F,u0,v0,pDOF)

%% Preliminaries
p = user_params(1);         % order of rational approximation
rhoInfty = user_params(2);  
N = p + 1;                  % no. of integration points
pf = N - 1;                 % order of force expansion

n = size(K,1);  % no. of degrees of freedom

K =  dt*dt*sparse(K);
C = dt*sparse(C);
M = sparse(M);
F =  dt*dt*F;
% Initial velocity --> normalized with dt
v0 = dt*v0;

%% Initialise output arrays
dsp = zeros(ns, length(pDOF)); % displacements
vel = dsp;  % velocities
acc = dsp;  % accelerations

%% Initialise scheme
[r, prcoe, cfr1] = InitSchemeRhoPF(p, rhoInfty, pf);

s = forceSamplingPoints(N);
Tcfr1 = transMtxPointsToPoly(s, p+1)*cfr1(1:p+1,:);

% Initial conditions
tm = 0;
z = [v0; u0]; 

% Effective stiffness
Kd = sparse((r*r)*M + r*C + K);     
dKd = decomposition(Kd);

% Initial acceleration
if rhoInfty == 0 || nnz(z) == 0
    an = acc(1,:)';  % initial acceleration not required.
else
    MLumped  = sum(M,2);
    ftmp = -(C*v0 + K*u0); 
    an = (F*signal(0)+ftmp)./MLumped;
end

% Store initial output
it = 1;
dsp(it,:) = u0(pDOF);
vel(it,:) = v0(pDOF);
acc(it,:) = an;

%% Time-stepping algorithm
for it = 2:ns
    
    ts = tm + dt*s;
    Fp = F.*reshape(signal(ts), 1, []);

    tm = tm + dt;
    
    x = zeros(2*n,1);
    for ip = 1:p
        g = x + prcoe(ip)*z;
        rfri = Fp*Tcfr1(:,ip); % r*{fri}
        x(1:n) = dKd\(r*(M*g(1:n)) - K*g(n+1:end) + rfri);
        x(n+1:end) = (x(1:n) + g(n+1:end))/r;
    end
    z = prcoe(p+1)*z + x;
    an = prcoe(end)*an + (r*x(1:n) - g(1:n));

    % store responses for output
    vel(it,:) = z(pDOF,1);
    dsp(it,:) = z(n+pDOF,1);
    acc(it,:) = an;

end

vel = vel/dt;
acc = acc/(dt*dt);


end
