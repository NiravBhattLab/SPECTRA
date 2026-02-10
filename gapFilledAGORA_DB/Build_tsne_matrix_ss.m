clear
load('ConsUmodel');
ConsUmodel.subSystems(6740) = ConsUmodel.subSystems{6740};
p=dir('./AGORA_new_db');
p = {p(3:end).name}';
p = setdiff(p,{'GapfilledModelsMILP_180s'});

% getting the count of all the subsystems
n_ss = numel(unique(ConsUmodel.subSystems));
u_ss = unique(ConsUmodel.subSystems);

% building the tsne matrix for all the subsystems
tsne_mat_all_ss = zeros(numel(p),n_ss);
% building the tsne matrix for newly added rxns subsystems
tsne_mat_new_ss = zeros(numel(p),n_ss);

for i=1:numel(p)
    % loading the new db model
    load(['./AGORA_new_db/',p{i}])

    % loading the existing AGORA model
    m = load(['./AGORA2/',p{i}]);
    m = m.model;
    new_rxns = setdiff(model.rxns,m.rxns);

    % getting the count of all the subsystems in the new db model
    for j=1:n_ss
        tsne_mat_all_ss(i,j) = sum(ismember(model.rxns,ConsUmodel.rxns(ismember(ConsUmodel.subSystems,u_ss{j}))));
        tsne_mat_new_ss(i,j) = sum(ismember(new_rxns,ConsUmodel.rxns(ismember(ConsUmodel.subSystems,u_ss{j}))));
    end
end

% removing the columns with all zeros
ids1 =all(tsne_mat_all_ss==0,1);
ids1(1) = 1; % this is an empty cell
tsne_mat_all_ss(:,ids1) = [];
ids2 = all(tsne_mat_new_ss==0,1);
ids2(1) = 1; % this is an empty cell
tsne_mat_new_ss(:,ids2) = [];

% adding the row names from the vector p and saving the matrix as csv files
tsne_mat_all_ss_tbl = array2table(tsne_mat_all_ss,'VariableNames',u_ss(~ids1),'RowNames',p);
tsne_mat_new_ss_tbl = array2table(tsne_mat_new_ss,'VariableNames',u_ss(~ids2),'RowNames',p);
writetable(tsne_mat_all_ss_tbl,'tsne_mat_all_ss.csv','WriteRowNames',true);
writetable(tsne_mat_new_ss_tbl,'tsne_mat_new_ss.csv','WriteRowNames',true);