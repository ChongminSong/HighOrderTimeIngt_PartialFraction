function [prcoe] = shiftPolyCoefficients(pcoe,r)

M = size(pcoe,2) - 1;           %order of polynomial
zc = zeros(size(pcoe,1),1);     %vector of zeros
prcoe = pcoe;
for ii = M:-1:1
    prcoe(:,ii:end) = [prcoe(:,ii)+r*prcoe(:,ii+1) r*prcoe(:,ii+2:end) zc] ...
                  - [zc prcoe(:,ii+1:end)];
end

end