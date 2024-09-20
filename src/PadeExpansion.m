function [pcoe, qcoe, r] = PadeExpansion(M,rhoInfty)

% mixed-order Padé expansion (L=M-1, M)

% M: order
% rhoInfty: a weight factor. It equals to the asymptote of the spectrum at high frequency 

L = M-1;
[p1, q1] = PadeCoeff(M, M);
[p2, q2] = PadeCoeff(M, L);
pcoe = rhoInfty*p1 + (1-rhoInfty)*[p2 zeros(M-L)];
qcoe = rhoInfty*q1 + (1-rhoInfty)*q2;

%% roots (factorization)
r = roots(fliplr(qcoe));
r(imag(r)<-1.d-6) = [];
[~, idx] = sort(imag(r));
r = r(idx);

end

function [p, q] = PadeCoeff(M,L)
    fc = @(x) factorial(x);
    ii=0:L; 
    p=(fc(M+L-ii))./(fc(ii).*fc(L-ii)); 
    ii=0:M; 
    q=(fc(M+L-ii)).*((-1).^ii)./(fc(ii).*fc(M-ii))*fc(M)/fc(L); 
end