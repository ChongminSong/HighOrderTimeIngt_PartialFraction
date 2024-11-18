function  [dsp,vel,acc] = TimeSolverPF(user_params,signal,ns,dt,K,M,C,F,u0,v0,pDOF)

%% Preliminaries
p = user_params(1);         % order of rational approximation
rhoInfty = user_params(2);  
N = p + 1;                  % no. of integration points
pf = N - 1;                 % order of force expansion

n = size(K,1);

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
[rho, plr, rs, a, cfr] = InitSchemePadePF(p,rhoInfty,pf);

s = forceSamplingPoints(N);
Tcfr = transMtxPointsToPoly(s, p+1)*cfr(1:p+1,:);

% Initial conditions
tm = 0;
z0 = [v0; u0];

% Effective stiffness 
nCmplx = floor(p/2);
nReal = mod(p,2);
if nReal > 0
    r = rs(1);
    dKd1 = decomposition(sparse((r*r)*M + r*C + K));
end
if nCmplx > 0
    L = cell(nCmplx,1);   U = L; LUp = L; LUq = L;
    for ic = 1:nCmplx
        r = rs(ic+nReal);
        [L{ic},U{ic},LUp{ic},LUq{ic}] = lu(sparse((r*r)*M + r*C + K),'vector'); % more efficient than using decomposition function
    end
end

% Initial acceleration
if rhoInfty == 0 || (nnz(z0) == 0 && nnz(signal(0)) == 0)
    a0 = acc(1,:)';  % initial acceleration not required.
else
    MLumped  = sum(M,2);
    ftmp = -(C*v0 + K*u0); 
    a0 = (F*signal(0)+ftmp)./MLumped;
end

% Store initial output
it = 1;
dsp(it,:) = u0(pDOF);
vel(it,:) = v0(pDOF);
acc(it,:) = a0;

%% Time-stepping algorithm
for it= 2:ns

    ts = tm + dt*s;
    Fp = F.*reshape(signal(ts), 1, []);

    tm = tm + dt;

    z = rho*z0;
    an = rho*a0;

    % real roots
    if nReal > 0
        r = real(rs(1));
        g = real(plr(1))*z0;
        fri = Fp*(real(Tcfr(:,1)));
        ftmp = dKd1\(r*(M*g(1:n)) - K*g(n+1:end) + r*fri);
        x = [ftmp ; (ftmp+g(n+1:end))/r];
        z = z + real(a(1))*x;
        an = an + real(a(1))*(r*x(1:n) - g(1:n));
    end
    
    % complex roots
    if nCmplx > 0
        for ic = 1:nCmplx
            r = rs(ic+nReal);
            cg = plr(ic+nReal)*z0;
            cfri = Fp*(Tcfr(:,ic+nReal));
            ctmp = plr(ic+nReal)*(r*(M*z0(1:n)) - K*z0(n+1:end)) + r*cfri;
            tmp(LUp{ic},:) = U{ic}\(L{ic}\ctmp(LUq{ic},:));
            x = [tmp; (tmp + cg(n+1:end))/r ];
            z = z + 2*real(a(ic+nReal)*x);
            an = an + 2*real(a(ic+nReal)*(r*tmp - cg(1:n)));
        end
    end
    z0 = z;
    a0 = an;

    % store responses for output
    vel(it,:) = z0(pDOF,1);
    dsp(it,:) = z0(n+pDOF,1);
    acc(it,:) = a0;
    
end

vel = vel/dt;
acc = acc/(dt*dt);

end