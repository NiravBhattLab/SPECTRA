clear
model = get_core_met_toy_model();
a = sprintcc(model,1e-4);
tol=1e-4;
consModel = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),a)));
coreMets = [2,3,4];
weights = zeros(numel(consModel.rxns),1);
weights(1:9)=1;
nSol=10;
altSolMethod = 'coreDirection';
[ConsModel,LPS] = sprintcore(consModel,[],tol,[],coreMets,weights,nSol,altSolMethod,[],[],[]);