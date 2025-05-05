function [MicComModel,removedCoreRxns] = BuildCommunityModels(Folder_path, abbr, Umodel, UmodelMap, media, lbs, ubs, ConsiderOtherTranRxn, Cores, bioCore, TransferCore, relAbun, tol, weights, probType,solveTime)
% USAGE:
%   [MicComModel,removedCoreRxns] = BuildCommunityModels(Folder_path, abbr, Umodel, UmodelMap, media, lbs, ubs, ConsiderOtherTranRxn, Cores, bioCore, TransferCore, relAbun, tol, weights)
%
% INPUTS:
%   Folder_path:          Matlab cell listing the paths to all the microbial
%                         model structures in .mat format. Lower (*lb) and upper
%                         (*ub) bounds can be provided for any reaction in the 
%                         microbial models other than default bounds. Reaction 
%                         ID (*rxns) and metabolite ID (*mets) should be in same
%                         format as in the corresponding Umodel. 
%                         Metabolite IDs(*mets) should include the compartment 
%                         info(Eg: glc_D[e], pyr[c])
%   abbr:                 Matlab cell listing model abbrevations. All rxns and mets 
%                         will have this prefix. Must be same order as in Folder_path
%   Umodel:               Matlab cell listing the COBRA model structure of Universal 
%                         models. If only one model is provided, then it will 
%                         considered as universal model for all the microbial models,
%                         else the variable 'UmodelMap' has to be provided. 
%                         All the exchange reactions ID should begin with 'EX_'
%                         The following fields are required in the universal model:
%                           * S    - `m x n` Stoichiometric matrix
%                           * b    - `m x 1` change in concentration with time
%                           * c    - `n x 1` Linear objective coefficients
%                           * lb   - `n x 1` Lower bounds on net flux
%                           * ub   - `n x 1` Upper bounds on net flux
%                           * mets - metabolite IDs
%                           * rxns - reaction IDs
%   UmodelMap:            Array of size equal to number of microbial models. The 
%                         values define the ids of the correspoding Universal 
%                         model that has to be used
%   media:                matlab structure with fields
%                           *exc_rxns: list of exchange reactions (should be the union 
%                                      set of all the exchange reactions in the 
%                                      universal models)
%                           *lb:       Lower bounds of corresponding reactions in exc_rxns
%                           *ub:       Upper bounds of corresponding reactions in exc_rxns
%
% OPTIONAL INPUTS:
%   lbs:                  vector of lower bounds to all microbial biomass reactions (Default: zeros) 
%   ubs:                  vector of upper bounds to all microbial biomass reactions (Default: 1000) 
%   ConsiderOtherTranRxn: Boolean value (Default: 1)
%                           1: All transport reactions that are present in the universal 
%                              model will be considered for inclusion
%                           0: Transport reactions that are present only in the microbial 
%                              model will only be considered for inclusion
%   Cores:                A cell of boolean vectors indicating the core reactions
%                         If Cores is provided TransferCore and bioCore
%                         will be ignored. (Default: {})
%   bioCore:              Boolen value (TransferCore will be ignored if bioCore is 1) (Default: 0)
%                           1: Only individual biomass reactions are considered as core
%                              reactions
%                           0: All the non-transport reactions in the individual models
%                              are considered as core reactions
%   TransferCore:         Boolen value (Default: 0)
%                           1: Transport reactions that are present in the microbial
%                              model will also be considered as core reactions
%                           0: Transport reactions that are present in the microbial
%                              model will not be considered as core reactions 
%   relAbun:              Vector denoting the relative abundance of all the strains
%                         order should be same as abbr (or) Folder_path
%                         (Default: ones)
%   tol:                  minimum absolute flux value that every reaction in the 
%                         community model should carry (Default: 1e-4)
%   weights:              A matlab cell of weight vectors. It should be same size as that
%                         of Umodel. More the weight for a reaction lesser
%                         the chances of getting included in the final
%                         model (Default: one to all the reactions)
%   probType:             Which optimization to use to find the minimal reaction set.
%                         accepted values: 'LP','MILP','DC'. (Default: LP)
%   solveTime:            Maximum runtime for solving MILP problem (Default: 7200s)
%
% OUTPUTS:
%   ComModel:             The consistent community model
%   removedRxns:          Reactions that were removed because of inconsistency 
%                         (or) inability to carry the minimum flux (tol) in the given
%                         community and media conditions
% 
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if numel(Folder_path) ~= numel(abbr)
    error('Model names has to be of the same size as the number of models in Foler_path')
end
if numel(UmodelMap) ~= numel(Folder_path)
    error('UmodelMap has to be of the same size as the number of models in Foler_path')
end

if ~exist('lbs', 'var') || isempty(lbs)
    lbs = zeros(numel(abbr),1);
end
if ~exist('ubs', 'var') || isempty(ubs)
    ubs = ones(numel(abbr),1)*1000;
end
if ~exist('ConsiderOtherTranRxn', 'var') || isempty(ConsiderOtherTranRxn)
    ConsiderOtherTranRxn = 1;
end
if ~exist('Cores', 'var') || isempty(Cores)
    Cores = {};
end
if ~exist('bioCore', 'var') || isempty(bioCore)
    bioCore = 0;
end
if ~exist('TransferCore', 'var') || isempty(TransferCore)
    TransferCore = 0;
end
if ~exist('relAbun', 'var') || isempty(relAbun)
    relAbun = ones(numel(abbr),1);
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-4;
end
if ~exist('weights', 'var') || isempty(weights)
    weights={};
    for i=1:numel(Umodel)
        temp = Umodel{i};
        weights{i,1} = ones(numel(temp.rxns),1);
    end
end
if ~exist('probType', 'var') || isempty(probType)
    probType='LP';  
end

n_models = numel(Folder_path); % number of models
n_Umodels = numel(Umodel); % number of universal models
Umodel_new = Umodel;

exc_rxns={};exc_rxnFormulas={};
for j=1:n_Umodels
    Utemp = Umodel_new{j};
    Curr_exc=Utemp.rxns(startsWith(Utemp.rxns,'EX_'));
    exc_rxns=[exc_rxns;Curr_exc];
    form =printRxnFormula(Utemp,'rxnAbbrList',Curr_exc,'printFlag',false);
    exc_rxnFormulas = [exc_rxnFormulas;form];
    Utemp=removeRxns(Utemp,exc_rxns); % remove all exchange rxns
    Umodel_new{j} = Utemp;
end

[~,ia,~] = unique(exc_rxns);
exc_rxns = exc_rxns(ia);
exc_rxnFormulas = exc_rxnFormulas(ia);

S=[];lb=[];ub=[];c=[];b=[];rxns=[];mets=[];core=[];weight=[];metNames=[];rxnNames=[];
for i=1:n_models
    load(Folder_path{i})
    CuUmodel = Umodel_new{UmodelMap(i)};
    if ConsiderOtherTranRxn
        UmodelTemp = CuUmodel;
    else
        % removing all the transport reactions that are not there in the
        % given microbial model
        UmodelTemp = removeRxns(CuUmodel,setdiff(CuUmodel.rxns(getTransRxns(CuUmodel)),...
            model.rxns(getTransRxns(model))));
    end
    
    if ~isempty(Cores)
        currCore = Cores{i};
        if numel(currCore)~= numel(model.rxns)
            error('Boolean vector indicating the core reactions must be of same size as the model')
        end
        coreTemp = zeros(numel(UmodelTemp.rxns),1);
        [ia,ib] = ismember(UmodelTemp.rxns,model.rxns);
        coreTemp(ia) = currCore(ib(ib~=0));
    elseif bioCore
        coreTemp = zeros(numel(UmodelTemp.rxns),1); % all reactions are non-core except biomass reactions
    elseif TransferCore
        coreTemp = ismember(UmodelTemp.rxns,model.rxns);
    else
        coreTemp = ismember(UmodelTemp.rxns,model.rxns(setdiff([1:numel(model.rxns)],getTransRxns(model))));
    end
    
    % Adding the biomass reaction
    BioForm = printRxnFormula(model,model.rxns(find(model.c)),0);
    if ~isempty(BioForm)
        UmodelTemp=addReaction(UmodelTemp,'bio','reactionFormula',BioForm{1},...
            'lowerBound',lbs(i),'upperBound',ubs(i));
        core = [core;coreTemp;1]; % one is for the biomass reaction
    else
        core = [core;coreTemp];
    end
    % Adding relative abundance constraints
    trIDS = getTransRxns(UmodelTemp);
    UmodelTemp.lb(trIDS) = UmodelTemp.lb(trIDS)*relAbun(i);
    UmodelTemp.ub(trIDS) = UmodelTemp.ub(trIDS)*relAbun(i);
    
    
    new_rxns = cellfun(@(x)rename_rxns(x,abbr{i}),UmodelTemp.rxns,'uni',false);
    rxns = [rxns;new_rxns];
    rxnNames =[rxnNames;UmodelTemp.rxnNames];
    new_mets=cellfun(@(x)rename_mets(x,abbr{i}),UmodelTemp.mets,'uni',false);
    mets=[mets;new_mets];
    metNames=[metNames;UmodelTemp.metNames];
    S = blkdiag(S,UmodelTemp.S);c=[c;UmodelTemp.c];b=[b;UmodelTemp.b];
    new_lb = UmodelTemp.lb;new_ub = UmodelTemp.ub;
    [loca,locb] = ismember(model.rxns,UmodelTemp.rxns);
    locb = locb(locb~=0);
    new_lb(locb)=model.lb(loca);new_ub(locb)=model.ub(loca);
    lb=[lb;new_lb];ub=[ub;new_ub];
    
    temp = Umodel{UmodelMap(i)};
    wt = weights{UmodelMap(i)};
    wtTemp = ones(numel(UmodelTemp.rxns),1);
    [ia,ib] = ismember(UmodelTemp.rxns,temp.rxns);
    wtTemp(ia) = wt(ib(ib~=0));
    weight = [weight;wtTemp];
end

% merging the extracellular metabolite rows and removing the extra
% metabolites
[Umets,~,ix] = unique(mets);
counts = accumarray(ix,1).';
counts = counts';
for j=1:numel(counts)
    if counts(j)>1
        ids = find(ismember(mets,Umets{j}));
        S(ids(1),:)=sum(S(ids,:),1);
        S(ids(2:end),:)=[];
        b(ids(2:end))=[];
        mets(ids(2:end))=[];
        metNames(ids(2:end))=[];
    end
end
    
ComModel=struct();
ComModel.S=S;ComModel.lb=lb;
ComModel.ub=ub;ComModel.c=c;
ComModel.b=b;ComModel.mets=mets;ComModel.metNames=metNames;
ComModel.rxns=rxns; ComModel.rxnNames=rxnNames;

for i=1:numel(media.exc_rxns)
    id = find(ismember(exc_rxns,media.exc_rxns{i}));
    if ~isempty(id)
        ComModel=addReaction(ComModel,exc_rxns{id},'reactionFormula',exc_rxnFormulas{id},...
            'lowerBound',media.lb(i),'upperBound',media.ub(i));
    end
end
core = [core;zeros(sum(ismember(exc_rxns,media.exc_rxns)),1)]; % exchange reactions are non-core
weight = [weight;zeros(sum(ismember(exc_rxns,media.exc_rxns)),1)];

if strcmp(probType,'MILP')
    if ~exist('solveTime', 'var') || isempty(solveTime)
        solveTime=7200;     
    end
    [MicComModel,removedCoreRxns] = SprintGapFiller(ComModel,core,weight,tol,1,probType,solveTime);
else
    [MicComModel,removedCoreRxns] = SprintGapFiller(ComModel,core,weight,tol,1,probType);
end



end


function rxns = rename_rxns(a,ABR)
rxns = [ABR,'_',a];
end

function mets = rename_mets(a,ABR)
if ~strcmp(a(end-1),'e')&&~contains(a,'biomass')
    mets = [ABR,'_',a];
else
    mets = a;
end
end

function ids = getTransRxns(model)
% this part of code is adapted from findTransRxns.m in COBRA toolbox
[~,compSymbols]=arrayfun(@parseMetNames, model.mets);
ids=[];
for i = 1:numel(model.rxns)
    % get compartment symbols for each rxn
    compSymbolsTmp=compSymbols(model.S(:,i)~=0);
    % if there's more than 1 compartment involved, it's a transport rxn
    if length(unique(compSymbolsTmp))>1
        ids=[ids;i];
    end
end
end
