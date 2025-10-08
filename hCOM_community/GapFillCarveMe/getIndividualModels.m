function models = getIndividualModels(Cmodel,modelNames)
% USAGE:
%   models = getIndividualModels(Cmodel,modelNames)
%
% INPUTS:
%    Cmodel:     The community model
%    modelNames: Names of the indivudual models. All the reactions and
%                metaboites in Cmodel has to be named after modelNames
%
% OUTPUTS:
%    models:     A cell consisting of the individual models in the order
%                same as that of modelNames

models={};
% adding the model specific reactions
for k =1:numel(modelNames)
    rxn_ids = find(startsWith(Cmodel.rxns,modelNames{k}));
    m = removeRxns(Cmodel,Cmodel.rxns(setdiff([1:numel(Cmodel.rxns)],rxn_ids)));
    models{k,1} = m;
end
% adding the exchange reactions based on the presence of the transport reactions
for j = 1:numel(modelNames)
    m = models{j};
    [~,compSymbols]=arrayfun(@parseMetNames, m.mets);
    id_ex_met =find(ismember(compSymbols,'e'));
    m = addExchangeRxn(m,m.mets(id_ex_met)); % adding exchange reactions for all the extracellular metabolites
    m = addExchangeRxn(m,'biomass[c]'); % adding exchange reaction for biomass (this is intracellular)
    m.rxns = strrep(m.rxns,[modelNames{j},'_'],''); % removing the microbe name in rxns
    m.mets = strrep(m.mets,[modelNames{j},'_'],''); % removing the microbe name in mets
    id = find(ismember(m.rxns,'bio')); % getting the objective function
    if numel(id)~=1
        warning(['No biomass reaction for ',modelNames{j}])
    else
        m.c(id)=1;
    end
    
    m.lb(m.lb<0)=-1000;
    m.ub(:)=1000;
    m.modelName= modelNames{j};
    models{j,1} = m;
end
end