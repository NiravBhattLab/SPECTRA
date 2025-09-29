clear
changeCobraSolverParams('LP', 'feasTol', 1e-9);
tbl = readtable('AGORA2_details.xlsx');
tbl2 = tbl(:,[1,13]);
nRxnsAG2=[];nRxnsAG3=[];
nBlkAG2=[];nBlkAG3=[];
p=dir('./AGORA2');
p={p(3:end).name}';
for i =1:numel(p)
    load(['./AGORA2/',p{i}])
    nRxnsAG2(i) = numel(model.rxns);
    a = spectraCC(model,1e-8);
    nBlkAG2(i) = numel(model.rxns)-numel(a);
    m = model;
    clear model
    load(['./AGORA3/GapfilledModelsMILP_180s/',p{i}])
    temp = numel(setdiff(m.rxns,model.rxns)); % reactions that could not be made consistent
    nRxnsAG3(i) = numel(model.rxns)+temp;
    nBlkAG3(i) = temp;
    clear model
end

ids=[];
asd= strrep(p,'.mat','');
for i=1:numel(p)
    ids(i) = find(ismember(asd,tbl2.MicrobeID{i}));
end

tbl2.AGORA2_nRxns = nRxnsAG2(ids)';
tbl2.AGORA3_nRxns = nRxnsAG3(ids)';
tbl2.nBlkdAG2 = nBlkAG2(ids)';
tbl2.nBlkdAG3 = nBlkAG3(ids)';
temp ={};
for i=1:7302
    temp{i,1} = [tbl2.MicrobeID{i},'|',num2str(tbl2.NCBITaxonomyID(i))];
end
tbl2.NodeID = temp;
save('Results_CountOfBlkdInAG2_AG3')