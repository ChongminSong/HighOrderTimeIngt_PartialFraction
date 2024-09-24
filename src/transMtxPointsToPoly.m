function [T] = transMtxPointsToPoly(s, nf)

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