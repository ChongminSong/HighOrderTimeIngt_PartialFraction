function [T] = transMtxPointsToPoly(s, nf)

% s = sample points
% nf = pf + 1

s  = reshape(s,[],1);
if length(s) >= nf
    a = (s-0.5).^(0:nf-1); % V^T
    T = a/(a'*a);
else
    disp([' ******* The number of sampling points: ', num2str(length(s))]);
    disp(['         should be larger than the number of terms of polynomial: ' ...
        num2str(nf)]);
end