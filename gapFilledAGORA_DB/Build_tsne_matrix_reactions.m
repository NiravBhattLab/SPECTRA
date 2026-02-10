clear
load('ConsUmodel');
p=dir('./AGORA_new_db');
p = {p(3:end).name}';
p = setdiff(p,{'GapfilledModelsMILP_180s'});
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
ids1 =all(tsne_mat_all_rxns==0,1);
tsne_mat_all_rxns(:,ids1) = [];
ids2 = all(tsne_mat_new_rxns==0,1);
tsne_mat_new_rxns(:,ids2) = [];

% adding the row names from the vector p and saving the matrix as csv files
tsne_mat_all_rxns_tbl = array2table(tsne_mat_all_rxns,'VariableNames',ConsUmodel.rxns(~ids1),'RowNames',p);
tsne_mat_new_rxns_tbl = array2table(tsne_mat_new_rxns,'VariableNames',ConsUmodel.rxns(~ids2),'RowNames',p);
writetable(tsne_mat_all_rxns_tbl,'tsne_mat_all_rxns.csv','WriteRowNames',true);
writetable(tsne_mat_new_rxns_tbl,'tsne_mat_new_rxns.csv','WriteRowNames',true);
