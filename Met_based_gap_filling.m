clear
model = get_core_met_toy_model();
weights = zeros(numel(model.rxns),1);
weights(1:9)=1;
nSol=10;
altSolMethod = 'pathwayExclusion';
coreRxns=[24,25,26];
tol =1e-4;
gapFilltype='topology';
probType ='MILP';
[Model,BlockedCoreRxns,LPs] = SprintGapFiller(model,coreRxns,tol,gapFilltype,weights,nSol,altSolMethod,probType);