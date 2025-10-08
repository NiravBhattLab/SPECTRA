clear
changeCobraSolver('gurobi','all')
changeCobraSolverParams('LP', 'feasTol', 1e-8);
changeCobraSolverParams('LP', 'optTol', 1e-8);
tol=1e-5;
load('Bigg2VMHEx.mat') % to rename the scores obtained

Folder_path = dir('./CarvemeDraftModels');
Folder_path = {Folder_path(3:end).name};
Folder_path = strrep(Folder_path,'.mat','');
abbr = Folder_path; % abbr defines the model names that will be added as prefix to the reations and metabolites

% Loading the universal models
load('./Universal_models/carvemeModelConv') % Carveme's universal model
u = model;
load('./Universal_models/carvemeModelGNconv') % Carveme's universal gram negative model
load('./Universal_models/carvemeModelGPconv') % Carveme's universal gram positive model

media = getSAACMedia(u,GrPosModel,GrNegModel); % getting the SAAC media
stain_info = readtable('gram_stain.csv'); % getting the gram staining details

% paramters for spectraME
consType = 'stoichiometry';
nSol = 1;
altSolMethod = {};
probType = 'minNetMILP';

% defining the lower and upper bounds for biomass reactions
lbs = 0.1*ones(numel(abbr),1);
ubs = 1000*ones(numel(abbr),1);

Umodel = {u,GrPosModel,GrNegModel}; % universal models in cell

n_models = numel(Folder_path); % number of models
n_Umodels = numel(Umodel); % number of universal models
Umodel_new = Umodel;

% removing all the exchange reactions from the universal models and storing them in the variable 'Umodel_new'
% Also storing the exchange reactions and their formulas in the variables 'exc_rxns' and 'exc_rxnFormulas'
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

% creating the universal model for the community
S=[];lb=[];ub=[];c=[];b=[];rxns=[];mets=[];core=[];weight=[];metNames=[];rxnNames=[];
for i=1:n_models
    load(Folder_path{i}) % loading the microbial model

    % setting the lower bound of ATPM to 0.1
    model.lb(ismember(model.rxns,'ATPM')) = 0.1;

    % mapping the universal model based on the gram staining
    stain_id = find(ismember(stain_info.Strain,abbr{i}));
    if strcmp(stain_info.stain{stain_id},'0')
        CuUmodel = Umodel_new{1};
    elseif strcmp(stain_info.stain{stain_id},'+')
        CuUmodel = Umodel_new{2};
    elseif strcmp(stain_info.stain{stain_id},'-')
        CuUmodel = Umodel_new{3};
    else
        error('Gram staining information not available');
    end

    UmodelTemp = CuUmodel;
    coreTemp = zeros(numel(UmodelTemp.rxns),1); % all reactions are non-core except biomass reactions
    
    % Adding the biomass reaction
    BioForm = printRxnFormula(model,model.rxns(find(model.c)),0);
    if ~isempty(BioForm)
        UmodelTemp=addReaction(UmodelTemp,'bio','reactionFormula',BioForm{1},...
            'lowerBound',lbs(i),'upperBound',ubs(i));
        UmodelTemp.c(end)=1; % setting the biomass reaction as objective
        core = [core;coreTemp;1]; % one is for the biomass reaction
    else
        core = [core;coreTemp];
    end

    % Adding relative abundance constraints
    % trIDS = getTransRxns(UmodelTemp);
    % UmodelTemp.lb(trIDS) = UmodelTemp.lb(trIDS)*relAbun(i);
    % UmodelTemp.ub(trIDS) = UmodelTemp.ub(trIDS)*relAbun(i);
    
    new_rxns = cellfun(@(x)rename_rxns(x,abbr{i}),UmodelTemp.rxns,'uni',false); % all the reactions will have the model abbr as prefix
    rxns = [rxns;new_rxns];
    rxnNames =[rxnNames;UmodelTemp.rxnNames];
    new_mets=cellfun(@(x)rename_mets(x,abbr{i}),UmodelTemp.mets,'uni',false); % all the metabolites except extracellular and biomass will have the model abbr as prefix
    mets=[mets;new_mets];

    metNames=[metNames;UmodelTemp.metNames];
    S = blkdiag(S,UmodelTemp.S);c=[c;UmodelTemp.c];b=[b;UmodelTemp.b]; % concatenating the S,c and b matrices
    
    % the bounds will be same as that in the microbial model
    % except for the reactions that are not present in the microbial model. For these reactions universal model constraints will be used
    new_lb = UmodelTemp.lb;new_ub = UmodelTemp.ub; 
    [loca,locb] = ismember(model.rxns,UmodelTemp.rxns);
    locb = locb(locb~=0);
    new_lb(locb)=model.lb(loca);new_ub(locb)=model.ub(loca);
    lb=[lb;new_lb];ub=[ub;new_ub];

    % getting the weights for the gapfilling step
    scores = readtable(['Carveme_results\',Folder_path{i},'\scores.tsv'],'FileType', 'text', 'Delimiter', '\t');
    scores = rename_reactions(scores,Bigg2VMHEx); % renaming the reactions in accordance with the universal model
    weights = ones(numel(UmodelTemp.rxns),1);
    temp_wts = 1./(1+scores.normalized_score);
    [ia,ib] = ismember(scores.reaction, UmodelTemp.rxns);
    weights(ib(ia)) = temp_wts(ia);
    % assigining zero weights to the existing reactions
    idx = ismember(UmodelTemp.rxns,model.rxns);
    weights(idx) = 0;
    weight = [weight;weights];
end


% merging the extracellular metabolite rows and removing the extra metabolites
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

% the universal model of the microbial community
ComModel=struct();
ComModel.S=S;ComModel.lb=lb;
ComModel.ub=ub;ComModel.c=c;
ComModel.b=b;ComModel.mets=mets;ComModel.metNames=metNames;
ComModel.rxns=rxns; ComModel.rxnNames=rxnNames;

% defining the SAAC media conditions. first we will collect all the exchange reactions present in all the universal models, then apply saac media constraints.
media=struct();
ids =startsWith(u.rxns,'EX_');
media.exc_rxns = u.rxns(ids);
media.lb = u.lb(ids); media.ub = u.ub(ids);
ids = startsWith(GrPosModel.rxns,'EX_');
media.exc_rxns = [media.exc_rxns;GrPosModel.rxns(ids)];
media.lb = [media.lb;GrPosModel.lb(ids)];
media.ub = [media.ub;GrPosModel.ub(ids)];
ids = startsWith(GrNegModel.rxns,'EX_');
media.exc_rxns = [media.exc_rxns;GrNegModel.rxns(ids)];
media.lb = [media.lb;GrNegModel.lb(ids)];
media.ub = [media.ub;GrNegModel.ub(ids)];
[C,ia,~] = unique(media.exc_rxns);
media.exc_rxns=media.exc_rxns(ia);
media.ub=media.ub(ia);
media.lb=media.lb(ia);

% SAAC media constraints
tbl = readtable('../SAACmedia.xlsx','Sheet','carveme_bounds');
tbl.ExchangeReactions = strrep(tbl.ExchangeReactions,'''','');
media.lb(:)=0;
for i =1:numel(tbl.ExchangeReactions)
    r=tbl.ExchangeReactions{i};
    b=tbl.Bounds(i);
    if b*100<10
        media.lb(ismember(media.exc_rxns,r))=-10;
    else
        media.lb(ismember(media.exc_rxns,r))=-b*100;
    end
end
media.ub=media.ub*100;% scaling the upper bound thereby allowing large-scale secretion


% appending the exchange reactions
for i=1:numel(media.exc_rxns)
    id = find(ismember(exc_rxns,media.exc_rxns{i}));
    if ~isempty(id)
        ComModel=addReaction(ComModel,exc_rxns{id},'reactionFormula',exc_rxnFormulas{id},...
            'lowerBound',media.lb(i),'upperBound',media.ub(i));
    end
end
core = [core;zeros(sum(ismember(exc_rxns,media.exc_rxns)),1)]; % exchange reactions are non-core
weight = [weight;ones(sum(ismember(exc_rxns,media.exc_rxns)),1)]; % providing the weights as ones to the exchange reactions

% consistency check to the universal community model
consistent_ids = spectraCC(ComModel,tol);
ConsComModel = removeRxns(ComModel,ComModel.rxns(setdiff(1:numel(ComModel.rxns),consistent_ids)));  

% updating the cores and weights accordingly
[~,ib] = ismember(ConsComModel.rxns,ComModel.rxns);
weight = weight(ib); core = core(ib);
[MiComModel,LPS] = spectraME(ConsComModel,find(core),tol,consType,weight,nSol,altSolMethod,probType,7200*5); % gapfilling the community model
models = getIndividualModels(MiComModel,abbr) % getting the individual models from the community model

% merging the actual models with the newly obtained models
for k=1:numel(abbr)
    load(['./CarvemeDraftModels/',abbr{k},'.mat']) % loading the microbial model
    gapFilledModel = models{k};
    Model = getUnionModel(gapFilledModel,model); 
    save(['./CarvemeCommunityGapfilling/',abbr{k}],'Model')
    clear model Model gapFilledModel
end

end

function Model = getUnionModel(Model,model)
    Model = mergeTwoModels(Model, model);
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
function scores = rename_reactions(scores,Bigg2VMHEx)
    [ia,ib] = ismember(scores.reaction,Bigg2VMHEx.rxns.bigg);
    scores.reaction(ia) = Bigg2VMHEx.rxns.vmh(ib(find(ib)));
    % conversion of carveme rxns
    scores.reaction = regexprep(scores.reaction,'^R_','');
    scores.reaction = regexprep(scores.reaction,'_e$','(e)');
    scores.reaction = regexprep(scores.reaction,'_c$','(c)');
end