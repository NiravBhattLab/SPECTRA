clear
load('./UpdatedRecon3D.mat')
tbl = readtable('../../recon-store-genes-1.tsv','FileType','text');
t = {};
for i=1:numel(tbl.gene_number)
    t{i,1} = num2str(tbl.gene_number(i));
end
new_genes = {};
for i =1:numel(model.genes)
    id = find(ismember(t,model.genes{i}));
    if isempty(id)
       new_genes{i,1} = model.genes{i};
    else
        new_genes{i,1} = tbl.symbol{i};
    end
end
model.genes = new_genes;
save("consRecon3DGeneSymbol",'model')