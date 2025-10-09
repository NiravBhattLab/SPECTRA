% gapfilling the individual hcom models using SAAC media
clear
changeCobraSolver('gurobi','all')
changeCobraSolverParams('LP', 'feasTol', 1e-8);
changeCobraSolverParams('LP', 'optTol', 1e-8);
tol=1e-5;
load('Bigg2VMHEx.mat')
p = dir('./CarvemeDraftModels/');
p = {p(3:end).name}';
p=  strrep(p,'.mat','');
% Loading the universal models
load('./Universal_models/carvemeModelConv') % Carveme's universal model
u = model;
load('./Universal_models/carvemeModelGNconv') % Carveme's universal gram negative model
load('./Universal_models/carvemeModelGPconv') % Carveme's universal gram positive model

media = getSAACMedia(u,GrPosModel,GrNegModel); % getting the SAAC media
stain_info = readtable('gram_stain.csv'); % getting the gram staining details

% paramters for spectraMe
consType = 'stoichiometry';
nSol = 1;
altSolMethod = {};
probType = 'minNetMILP';

p2 = dir('./CarvemeIndividualGapfilling/');
p2 = {p2(3:end).name}';
p2 = strrep(p2,'.mat','');

remaining = setdiff(p,p2);

for k=1:numel(remaining) 
    try
        i = find(ismember(stain_info.Strain,remaining{k}));
        load(['./CarvemeDraftModels/',p{i}])
        if strcmp(stain_info.stain{i},'0')
            Umodel1 = u;
        elseif strcmp(stain_info.stain{i},'+')
            Umodel1 = GrPosModel;
        elseif strcmp(stain_info.stain{i},'-')
            Umodel1 = GrNegModel;
        else
            error('Gram staining information not available');
        end
        ids =startsWith(Umodel1.rxns,'EX_'); % list of exchange reactions
        [ia,ib] = ismember(Umodel1.rxns(ids),media.exc_rxns);
        if sum(ia)~=numel(ia)
            error('Exchange reactions in the universal model not found in the media');
        end
        Umodel1.lb(ids) = media.lb(ib(ia));
        Umodel1.ub(ids) = media.ub(ib(ia));
        % adding the biomass reaction to the universal model
        bio_form = printRxnFormula(model,model.rxns(find(model.c)),0);
        bio_rxn = model.rxns(find(model.c));
        bio_rxnName = model.rxnNames(find(model.c));

        Umodel1 = addReaction(Umodel1,bio_rxn{1},'reactionName',bio_rxnName{1},'reactionFormula',bio_form{1});
        Umodel1.c(:)=0;
        Umodel1.c(ismember(Umodel1.rxns,bio_rxn{1}))=1;

        % building a consistent universal model
        ConsReacIDS = spectraCC(Umodel1,tol);
        if ~ismember(find(Umodel1.c),ConsReacIDS)
            error('inconsistent biomass reaction')
        end
        Umodel = removeRxns(Umodel1,Umodel1.rxns(setdiff(1:numel(Umodel1.rxns),ConsReacIDS)));
        scores = readtable(['Carveme_results\',p{i},'\scores.tsv'],'FileType', 'text', 'Delimiter', '\t');
        scores = rename_reactions(scores,Bigg2VMHEx); % renaming the reactions in accordance with the universal model

        % verification step
        if sum(ismember(scores.reaction,Umodel1.rxns))~=numel(scores.reaction)
            error('missing reactions in the universal model')
        end

        % getting the weights for the gapfilling step
        weights = ones(numel(Umodel.rxns),1);
        temp_wts = 1./(1+scores.normalized_score);
        [ia,ib] = ismember(scores.reaction, Umodel.rxns);
        weights(ib(ia)) = temp_wts(ia);
        % assigining zero weights to the existing reactions
        idx = ismember(Umodel.rxns,model.rxns);
        weights(idx) = 0;

        % lower bound to the biomas reaction
        Umodel.lb(ismember(Umodel.rxns,bio_rxn{1})) = 0.1;
        % lower bound to the ATPM reaction
        Umodel.lb(ismember(Umodel.rxns,'ATPM')) = 0.1;
        [Model,LPS] = spectraME(Umodel,[],tol,consType,weights,nSol,altSolMethod,probType,7200);
        Model = getUnionModel(Model,model);
        % chekching if the model grows
        sol = optimizeCbModel(Model);
        if sol.f ==0 
            error('The microbe does not grow in the provided media')
        end
        save(['./CarvemeIndividualGapfilling/',p{i}],'Model')

        clear model Umodel scores weights
    catch ME
        disp(['Error in processing the model ',remaining{k}])
        disp(ME.message)
        continue
    end
end


function media = getSAACMedia(u,GrPosModel,GrNegModel)
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
    
    [~,ia,~] = unique(media.exc_rxns);
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
end

function scores = rename_reactions(scores,Bigg2VMHEx)
    [ia,ib] = ismember(scores.reaction,Bigg2VMHEx.rxns.bigg);
    scores.reaction(ia) = Bigg2VMHEx.rxns.vmh(ib(find(ib)));
    % conversion of carveme rxns
    scores.reaction = regexprep(scores.reaction,'^R_','');
    scores.reaction = regexprep(scores.reaction,'_e$','(e)');
    scores.reaction = regexprep(scores.reaction,'_c$','(c)');
end

function Model = getUnionModel(Model,model)
    Model = mergeTwoModels(Model, model);
end