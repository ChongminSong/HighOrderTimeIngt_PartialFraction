function [rho, plr, rs, a, cfr] = InitSchemePadePF(M, rhoInfty, pf)

[pcoe, qcoe, rs] = PadeExpansion(M,rhoInfty);

rho = pcoe(end)/qcoe(end);
a = polyPartialFraction(qcoe, rs);

plcoe = pcoe(1:end-1) - rho*qcoe(1:end-1);
plr = plcoe*reshape(rs,1,[]).^((0:M-1)');

cf = TimeIntgCoeffForce(pcoe,qcoe,pf);
cfr = cf*reshape(rs,1,[]).^((0:M-1)');

end