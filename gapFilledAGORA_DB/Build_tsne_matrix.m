clear
load('ConsUmodel');
p=dir('./AGORA_new_db');
p = {p(3:end).name}';
% building the tsne matrix for all the rxns
tsne_mat_all_rxns = zeros(numel(p),numel(ConsUmodel.rxns));
% building the tsne matrix for newly added rxns
tsne_mat_new_rxns = zeros(numel(p),numel(ConsUmodel.rxns));
for i=1:numel(p)
    % loading the new db model
    load(['./AGORA_new_db/',p{i}])
    tsne_mat_all_rxns(i,:) = ismember(ConsUmodel.rxns,model.rxns);

    % loading the existing AGORA model
    m = load(['./AGORA2/',p{i}]);
    m = m.model;
    new_rxns = setdiff(model.rxns,m.rxns);
    tsne_mat_new_rxns(i,:) = ismember(ConsUmodel.rxns,new_rxns);
end

% removing the columns with all zeros
tsne_mat_all_rxns(:,all(tsne_mat_all_rxns==0,1)) = [];
tsne_mat_new_rxns(:,all(tsne_mat_new_rxns==0,1)) = [];

% saving the matrices as csv files
writematrix(tsne_mat_all_rxns,'tsne_mat_all_rxns.csv')
writematrix(tsne_mat_new_rxns,'tsne_mat_new_rxns.csv')