clear
p = dir('./AGORA2');
p = {p(3:end).name}';

% to add the inconsistent reactions to the model
for i=1:numel(p)
    % agora db model
    m = load(['./AGORA2/',p{i}]);
    m =m.model;

    % new agora db model
    load(['./AGORA_new_db/GapfilledModelsMILP_180s/',p{i}]);
    model = removeRxns(model,intersect(m.rxns,model.rxns));
    % removing the unnecessary fields in m
    if isfield(m,'A')
        m = rmfield(m,'A');
    end
    if isfield(m,'C')
        m = rmfield(m,'C');
    end
    if isfield(m,'d')
        m = rmfield(m,'d');
    end
    model = mergeTwoModels(model,m);
    save(['./AGORA_new_db/',p{i}],'model')
end