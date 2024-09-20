function [s] = forceSamplingPoints(np)

% Gauss-Lobatto points
xi = lglnodes(np-1);
xi = flip(xi);

% Mapping to the interval [0,1]
scl = 1 - 1.0d-12;
s  = 1/2*(scl*xi + 1);

end