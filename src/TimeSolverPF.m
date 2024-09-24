function  [dsp,vel,acc] = TimeSolverPF(p,rhoInfty,signal,ns,dt,K,M,C,F,u0,v0,pDOF)

% The function is a part of the source code implementing 
% the high-order implicit time integration schemes in the paper:
% @misc{oshea2024highorderimplicittimeintegration,
%       title={A high-order implicit time integration method for linear and nonlinear dynamics with efficient computation of accelerations}, 
%       author={Daniel O'Shea and Xiaoran Zhang and Shayan Mohammadian and Chongmin Song},
%       year={2024},
%       eprint={2409.13397},
%       archivePrefix={arXiv},
%       primaryClass={math.NA},
%       url={https://arxiv.org/abs/2409.13397}, 
% }
% 
% The code is developed by:
%          Chongmin Song (c.song@unsw.edu.au)
%          Daniel O'Shea (d.oshea@unsw.edu.au)
%          Xiaoran Zhang (xiaoran.zhang1@student.unsw.edu.au)  
% 
% It is distributed under the MIT License (https://opensource.org/license/mit/)

%% Preliminaries
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
pf = p+1;
[rho, plr, rs, a, cfr] = InitSchemePadePF(p,rhoInfty,pf);

s = forceSamplingPoints(pf+1);
Tcfr = transMtxPointsToPoly(s, pf)*cfr(1:pf,:);

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
        y = [ftmp ; (ftmp+g(n+1:end))/r];
        z = z + real(a(1))*y;
        an = an + real(a(1))*(r*y(1:n) - g(1:n));
    end
    
    % complex roots
    if nCmplx > 0
        for ic = 1:nCmplx
            r = rs(ic+nReal);
            cg = plr(ic+nReal)*z0;
            cfri = Fp*(Tcfr(:,ic+nReal));
            ctmp = plr(ic+nReal)*(r*(M*z0(1:n)) - K*z0(n+1:end)) + r*cfri;
            tmp(LUp{ic},:) = U{ic}\(L{ic}\ctmp(LUq{ic},:));
            y = [tmp; (tmp + cg(n+1:end))/r ];
            z = z + 2*real(a(ic+nReal)*y);
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