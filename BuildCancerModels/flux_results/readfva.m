p=dir('../fva');
p = {p(1:end).name}';
p=p(find(contains(p,'ACH-')));

%% extract the flux variability results of each of cancer and non-cancerous models into a matlab structure
for i =1:length(p)
    fsrall=readtable(['../fva/',p{i,1}]);
    X=table2cell(fsrall);
    cells_all{:,i}=X;
end

save('cells_recon_2025',"cells_all")




