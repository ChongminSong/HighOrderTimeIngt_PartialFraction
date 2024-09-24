function [rho, plr, rs, a, cfr] = InitSchemePadePF(M, rhoInfty, pf)

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

[pcoe, qcoe, rs] = PadeExpansion(M,rhoInfty);

rho = pcoe(end)/qcoe(end);
a = polyPartialFraction(qcoe, rs);

plcoe = pcoe(1:end-1) - rho*qcoe(1:end-1);
plr = plcoe*reshape(rs,1,[]).^((0:M-1)');

cf = TimeIntgCoeffForce(pcoe,qcoe,pf);
cfr = cf*reshape(rs,1,[]).^((0:M-1)');

end