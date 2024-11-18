function [r] = MschemeRoot(M, rhoInfty)

RHS = [1,  1, -1,  1, -1, -1]*rhoInfty;
ir  = [1,  2,  2,  2,  3,  3];
j = 0:M;
pMcoe = ((-1).^j).*factorial(M)./factorial(j)./(factorial(M-j).^2);
pMcoe(end) = pMcoe(end) - RHS(M);
rs = sort( roots(pMcoe) );
r  = rs(ir(M));

end
