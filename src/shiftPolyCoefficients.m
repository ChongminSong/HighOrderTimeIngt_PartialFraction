function [prcoe] = shiftPolyCoefficients(pcoe,r)

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

M = size(pcoe,2) - 1;           %order of polynomial
zc = zeros(size(pcoe,1),1);     %vector of zeros
prcoe = pcoe;
for ii = M:-1:1
    prcoe(:,ii:end) = [prcoe(:,ii)+r*prcoe(:,ii+1) r*prcoe(:,ii+2:end) zc] ...
                  - [zc prcoe(:,ii+1:end)];
end

end