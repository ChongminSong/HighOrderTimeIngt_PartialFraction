function [dsp,vel,acc] = TimeSolverTRhoPF(p,rhoInfty,signal,ns,dt, ...
                                       K,M,C,F,u0,v0,pDOF)

%% Preliminaries
nz = size(K,1);

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
pf = p;
[r, prcoe, cfr1] = InitSchemeRhoPF(p, rhoInfty, pf);

s = forceSamplingPoints(pf+1);
Tcfr1 = transMtxPointsToPoly(s, pf+1)*cfr1(1:pf+1,:);

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
    
    zi = zeros(2*nz,1);
    for ip = 1:p
        g = zi + prcoe(ip)*z;
        rfri = Fp*Tcfr1(:,ip); % r*{fri}
        zi(1:nz) = dKd\(r*(M*g(1:nz)) - K*g(nz+1:end) + rfri);
        zi(nz+1:end) = (zi(1:nz) + g(nz+1:end))/r;
    end
    z = prcoe(p+1)*z + zi;
    an = prcoe(end)*an + (r*zi(1:nz) - g(1:nz));

    % store responses for output
    vel(it,:) = z(pDOF,1);
    dsp(it,:) = z(nz+pDOF,1);
    acc(it,:) = an;

end

vel = vel/dt;
acc = acc/(dt*dt);


end
