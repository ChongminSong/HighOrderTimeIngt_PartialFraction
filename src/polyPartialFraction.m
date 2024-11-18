function [a] = polyPartialFraction(q, r)

%% partial fraction of rational polynomial
M = length(q) - 1;
nterm = ceil(M/2);
nReal = mod(M,2);
a = zeros(nterm,1);
if M == 1
    a = 1;
else
    allRoots = [r;  conj(r(nReal+1:end))];
    for ii = 1:nterm
        a(ii) = 1/prod(allRoots([1:ii-1,ii+1:end])-r(ii));
    end
end

end