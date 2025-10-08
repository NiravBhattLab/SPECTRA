clear
p = dir("Carveme_results\");
p = {p(3:end).name}';
load('Bigg2VMHEx.mat')
for i=1:numel(p)
    model = readCbModel(['./Carveme_results/',p{i},'/model.xml']);
    model = rename_models_props(model,Bigg2VMHEx);
    save(['CarvemeDraftModels\',p{i}],'model')
    clear model
end

function model = rename_models_props(model,Bigg2VMHEx)
    [ia,ib] = ismember(model.mets,Bigg2VMHEx.mets.bigg);
    model.mets(ia) = Bigg2VMHEx.mets.vmh(ib(find(ib)));
    
    [ia,ib] = ismember(model.rxns,Bigg2VMHEx.rxns.bigg);
    model.rxns(ia) = Bigg2VMHEx.rxns.vmh(ib(find(ib)));
    
    % conversion of carveme mets
    model.mets = strrep(model.mets,'[C_c]','[c]');
    model.mets = strrep(model.mets,'[C_e]','[e]');
    model.mets = strrep(model.mets,'[C_p]','[p]');
    % conversion of carveme rxns
    model.rxns = regexprep(model.rxns,'^R_','');
    model.rxns = regexprep(model.rxns,'_e$','(e)');
    model.rxns = regexprep(model.rxns,'_c$','(c)');
    
    bio_form = printRxnFormula(model,model.rxns(find(model.c)),0);
    bio_rxn = model.rxns(find(model.c));
    bio_rxnName =model.rxnNames(find(model.c));
    model = removeRxns(model,model.rxns(find(model.c)));
    model = addReaction(model,bio_rxn{1},'reactionName',bio_rxnName{1},'reactionFormula',[bio_form{1},'+ biomass[c]']);
    model = addReaction(model,'EX_biomass(e)','reactionName','EX Biomass c0','reactionFormula','biomass[c]  <=> ');
    model.c(ismember(model.rxns,bio_rxn{1}))=1;
end