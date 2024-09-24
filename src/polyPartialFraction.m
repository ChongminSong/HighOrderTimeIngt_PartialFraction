function [a] = polyPartialFraction(q, r)

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