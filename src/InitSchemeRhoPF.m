function [r, prcoe, cfr] = InitSchemeRhoPF(M, rhoInfty, pf )

r = MschemeRoot(M, rhoInfty);

pcoe = pCoefficients(M, r);
prcoe = shiftPolyCoefficients(pcoe,r);

qrcoe = [zeros(1,M) 1];
qcoe = shiftPolyCoefficients(qrcoe,r);

cf = TimeIntgCoeffForce(pcoe,qcoe,pf);
cfr = r*shiftPolyCoefficients(cf,r);

end