% This code randomly deletes a % of reactions in AGORA models to get
% the incomplete models that needs to be gapfilled
clear
% initCobraToolbox(0);
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('LP', 'optTol', 1e-9);
AGORApath = './ConsAGORA';
p=dir(AGORApath);
p = {p(3:end).name}';
for i=1:numel(p)
    load([AGORApath,'/',p{i}])
    remPer = 0.3; % percentage of reactions to be removed
    n_rxns = numel(model.rxns);
    r = sort(randsample(n_rxns,int64(n_rxns*remPer)));
    r = setdiff(r,find(model.c)); % excluding the biomass reaction
    model = removeRxns(model,model.rxns(r));
    save(['./IncompModels_40/',p{i}],'model') 
end