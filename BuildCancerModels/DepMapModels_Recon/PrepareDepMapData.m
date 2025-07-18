clear
tbl =readtable('Expression_Public_23Q4.csv'); % data needs to be obtained from 
tbl = renamevars(tbl,'Var1','Samples');
geneExpression.value = table2array(tbl(:,2:end));
geneExpression.value = (2.^geneExpression.value)-1;
geneExpression.genes = tbl.Properties.VariableNames(2:end)';
geneExpression.context = tbl.Samples;
save('DepMapData','geneExpression')