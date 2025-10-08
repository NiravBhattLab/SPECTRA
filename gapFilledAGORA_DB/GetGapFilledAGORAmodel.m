% This code gap-fills the existing AGORA2 models using the
% available Universal AGORA model
clear
% initCobraToolbox(0);
changeCobraSolver('gurobi','all')
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'feasTol', 1e-9);
AGORApath = './AGORA2/';
p=dir(AGORApath);
p = {p(3:end).name}';
p2 = dir('./GapfilledModelsMILP_180s');
p2 = {p2(3:end).name}';
load('./ConsUmodel.mat')
p=setdiff(p,p2);
tol=1e-4;
for i=1:numel(p)
    load([AGORApath,'/',p{i}])
    n1=numel(model.rxns);
    % adding biomass reaction to the universal model
    bioR = printRxnFormula(model,'rxnAbbrList',model.rxns(find(model.c)),...
        'printFlag',false);
    Utemp = ConsUmodel;
    Utemp = addReaction(Utemp,model.rxns{find(model.c)},'reactionFormula',bioR{1},...
            'lowerBound',0,'upperBound',1000);
    % getting the core reactions
    core = ismember(Utemp.rxns,model.rxns);
    tic
    model = spectraCCME(Utemp,core,tol,'stoichiometry',[],1,{},'minNetMILP',180);
    t = toc;
    if t<60
        model = spectraCCME(Utemp,core,tol,'stoichiometry',[],1,{},'minNetMILP',180);
    end
    % checking if the model is consistent
    a = spectraCC(model,tol);

    if numel(a)==numel(model.rxns)
        % save the model if it is consistent
        save(['./GapfilledModelsMILP_180s/',p{i}],'model')
    else
        % save the model if inconsistent and also not it down
        save(['./GapfilledModelsMILP_180s/',p{i}],'model')
        load('inconsModels')
        inconsModels{end+1,1}=p{i};
        save('inconsModels','inconsModels')
        clear 'inconsModels'
    end
    % note on large sized models
    if numel(model.rxns)/n1>1.5
        load('HugeSizeModels')
        HugeSizeModels{end+1,1}=p{i};
        save('HugeSizeModels','HugeSizeModels')
        clear HugeSizeModels
    end
    
    clear model
end