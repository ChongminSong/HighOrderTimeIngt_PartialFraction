function [C] = TimeIntgCoeffForce(p, q, pf)

% coefficients for integration of non-homogeneous term
M = length(q) - 1;
tmp = p - q;
C = zeros(pf+1,M);
C(1,:) = tmp(2:end);
for k = 1:pf
    tmp = ((-1/2)^k)*(p-((-1)^k)*q);
    tmp(1:M) = tmp(1:M) + k*C(k,:); 
    C(k+1,:) =  tmp(2:end);
end