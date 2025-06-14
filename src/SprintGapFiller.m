function [Model,BlockedCoreRxns,flux,LPs] = SprintGapFiller(model,core,tol,gapFilltype,weights,nSol,altSolMethod,probType,solveTime,remGene)
% USAGE: 
%   [Model,BlockedCoreRxns,flux,LPs] = SprintGapFiller(model,core,tol,gapFilltype,weights,nSol,altSolMethod,probType,solveTime,remGene)
%
% INPUTS:
%     model:   COBRA model structure.
%     core:    core reactions which has to be present in the final model
%              (If any of the core reactions are blocked in the input model
%              then it will not be included in the final model and returned
%              as BlockedCoreRxns)
% 
% OPTIONAL INPUTS:
%     tol:          Minimum absolute flux required for a reaction to be unblocked (Default: 1e-4)
%     gapFilltype:  Type of gapfilling to apply. Either 'topology' or
%                   'stoichiometry' based. (Default:'stoichiometry')
%     weights:      Weights for non-core reactions. More the weights, lesser
%                   the chance to get included in the final model (Default: ones)
%     nSol:         Number of alternative solutions required (Default: 1)
%     altSolMethod: Method to find the alternate solutions.
%                   accepted values: 'coreDirection', 'pathwayExclusion'.
%                   Note: 'pathwayExclusion' works only for MILP probType
%                   (Default: 'coreDirection')
%     probType:     Which method to use to find the minimal reaction set.
%                   accepted values: 'LP','MILP','DC'. (Default: 'LP')
%     solveTime:    Maximum runtime for solving MILP problem (Default: 7200s)
%     remGene:      Bool value indicating whether to remove the unused genes
%                   or not (Default: 0 (doesn't remove the unused genes))
% OUTPUTS:
%     Model:           The consistent model (If nSol ==1). A cell
%                      consisting of models (If nSol>1).
%     BlockedCoreRxns: Core reactions that cannot have absolute flux 
%                      above the tol value in the input model
% 
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('weights', 'var') || isempty(weights)
    weights = ones(numel(model.rxns),1);    
end
if ~exist('tol', 'var') || isempty(tol)
    tol =1e-4;    
end
if ~exist('nSol', 'var') || isempty(nSol)
    nSol =1;
end
if ~exist('probType', 'var') || isempty(probType)
    probType='LP';  
end
if ~exist('remGene', 'var') || isempty(remGene)
    remGene=0;  
end
if ~exist('gapFilltype', 'var') || isempty(gapFilltype)
    gapFilltype='stoichiometry';  
end


[~,n] = size(model.S); % number of metabolites and reactions
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
rev = ones(n,1);
rev(model.lb>=0) = 0;
model.rev=rev;
prev_rxns = false(n,1);
temp_core = true(n,1);
flux = zeros(n,1); % initiating the flux vector that will carry directionality info
LPs=0;
core = ismember(1:n,core)';

if strcmp(gapFilltype,'stoichiometry')
    steadystate = 1;
elseif strcmp(gapFilltype,'topology')
    steadystate = 0;
end


while sum(temp_core)~=sum(prev_rxns)
    LPs = LPs+2;
    prev_rxns = temp_core;
    % maximizing number of reactions with forward flux
    [flux1,~] = forwardcc(model,temp_core,tol,steadystate);
    temp_core(abs(flux1)>=tol*0.99)=false;
    if sum(abs(flux))==0
        flux = flux1;
    else
        c1=round(unifrnd(0.45,0.55,1),4);
        flux = (c1*flux)+((1-c1)*flux1);
    end
    % maximizing number of reactions with reverse flux
    [flux2,~] = reverse(model,temp_core,tol,steadystate);
    temp_core(abs(flux2)>=tol*0.99)=false;
    
    c1=round(unifrnd(0.45,0.55,1),4);
    flux = (c1*flux)+((1-c1)*flux2);
    
end
BlckdRxns = find(temp_core);
BlockedCoreRxns = model.rxns(intersect(find(core),BlckdRxns));
core(temp_core) = 0; % not forcing any flux through the blocked reactions
direction = zeros(n,1);
direction(core==1&flux>0) = 1;
direction(core==1&flux<0) = -1;
if any(core==1&flux==0)
    % what if convex combination of a flux obtained at two iterations
    % cancel out each other
    warning('Any of the core reactions carry zero flux have to rerun again')
end

if strcmp(probType,'MILP')
    if ~exist('solveTime', 'var') || isempty(solveTime)
        solveTime=7200;     
    end
    [reacInd,x] = findConsistentReacID(model,direction,weights,tol,steadystate,probType,solveTime,flux);
elseif strcmp(probType,'LP')
    LPs=LPs+1;
    [reacInd,x] = findConsistentReacID(model,direction,weights,tol,steadystate,probType);
elseif strcmp(probType,'DC')
    [reacInd,x] = findConsistentReacID(model,direction,weights,tol,steadystate,probType);
end

flux = x(1:numel(model.rxns));
Model = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));
if remGene
    Model = removeUnusedGenes(Model);
end
if nSol>1
    % removing all the blocked reactions (this will be a consistent model)
    Nmodel = removeRxns(model,model.rxns(BlckdRxns));
    core2 = find(ismember(Nmodel.rxns,model.rxns(core)));
    [~,id] = ismember(Nmodel.rxns,model.rxns);
    weights2 = weights(id);
    reacInd2 = find(ismember(Nmodel.rxns,model.rxns(reacInd)));
    Models = sprintcore(Nmodel,core2,tol,weights2,nSol-1,altSolMethod,probType,solveTime,remGene,{reacInd2});
    Model=[Models;Model];
end

end