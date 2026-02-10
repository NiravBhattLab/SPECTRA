clear
model = get_core_met_toy_model();
weights = zeros(numel(model.rxns),1);
weights(1:9)=1;
nSol=10;
altSolMethod = 'pathwayExclusion';
coreRxns=[24,25,26];
tol =1e-4;
consType='topology';
probType ='minNetMILP';
[Model,BlockedCoreRxns,LPs] = spectraCCME(model,coreRxns,tol,consType,weights,nSol,altSolMethod,probType);