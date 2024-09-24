function [C] = TimeIntgCoeffForce(p, q, pf)

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