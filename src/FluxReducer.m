function flux = FluxReducer(model, steadyState)
% USAGE:
%     flux = FluxReducer(model)
%
% INPUTS:
%     model:       COBRA model structure.
%     steadyState: Boolean value indicating whether to assume steady state
%                  condition (S.v = 0) or accumulation condition (S.v >= 0)
%
% OUTPUTS:
%     flux:  A sparse flux vector obtained after reducing the sum of
%            absolute values of the fluxes
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

[m,n] = size(model.S);
rev = model.rev; % all the reversible reactions
n_ = sum(rev); % number of reversible reactions

% objective
f = [~rev;ones(n_,1)];

% equalities
Aeq = [model.S, sparse(m,n_)];
beq = zeros(m,1);
if steadyState
    csenseeq = repmat('E',m,1); % equality (For consistency based gap filling)
else
    csenseeq = repmat('G',m,1); % greater than (For topology based gap filling)
end


% inequalities
temp1 = speye(n);
temp2 = speye(n_);
Aineq1 = [temp1(rev,:),temp2];
bineq1 = zeros(n_,1);
csenseineq1 = repmat('G',n_,1); % greater than

Aineq2 = [temp1(rev,:),-1*temp2];
bineq2 = zeros(n_,1);
csenseineq2 = repmat('L',n_,1); % lesser than

% bounds
lb = model.lb;
lb = [lb;zeros(n_,1)];
ub = model.ub;
ub = [ub;Inf(n_,1)];

% Set up LP problem
LPproblem.A=[Aeq;Aineq1;Aineq2];
LPproblem.b=[beq;bineq1;bineq2];
LPproblem.lb=lb;
LPproblem.ub=ub;
LPproblem.c=f;
LPproblem.osense=1;%minimise
LPproblem.csense = [csenseeq; csenseineq1; csenseineq2];
solution = solveCobraLP(LPproblem);
stat =solution.stat;
if stat~=1
    fprintf('%s%s\n',num2str(solution.stat),' = solution.stat')
    fprintf('%s%s\n',num2str(solution.origStat),' = solution.origStat')
    warning('LP solution may not be optimal')
    flux=[];
else
    x=solution.full;
    flux = x(1:n);
end