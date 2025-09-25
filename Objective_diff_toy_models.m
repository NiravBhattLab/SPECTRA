clear
model = three_pathway_toy_model();
model.ub(:) =10;
model.ub([1,6,9]) = 1;

model.c(5)=1;
core = 5;
tol=1e-5;

% the LP based objective
consType = 'stoichiometry';
weights = ones(numel(model.rxns),1);
probType = 'minNetLP';
altSolMethod = {};
nSol=1;
[Model_minNetLP] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType);

% the minMILP based objective
consType = 'stoichiometry';
weights = ones(numel(model.rxns),1);
probType = 'minNetMILP';
altSolMethod = {};
nSol=1;
[Model_minNetMILP] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType);

% the DC based objective
consType = 'stoichiometry';
weights = ones(numel(model.rxns),1);
probType = 'minNetDC';
altSolMethod = {};
nSol=1;
[Model_minNetDC] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType);

% the tradeOff based objective
consType = 'stoichiometry';
weights = [1,1,-2,-1,1,1,2,-1,-2,0,1]';
probType = 'tradeOff';
altSolMethod = {};
nSol = 1;
[Model_tradeOff] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType);

% the optimBiomass based objective
consType = 'stoichiometry';
weights = ones(numel(model.rxns),1);
weights(find(model.c)) =10;
probType = 'growthOptim';
altSolMethod = {};
nSol = 1;
[Model_growthOptim] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType);

% to get the alternate solutions in MILP objectives

% the minMILP based objective
consType = 'stoichiometry';
weights = ones(numel(model.rxns),1);
probType = 'minNetMILP';
altSolMethod = 'pathwayExclusion';
nSol=5;
[Model_minNetMILP_altSols] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType);

% the tradeoff based objective
consType = 'stoichiometry';
weights = [1,1,-2,-1,1,1,2,-1,-2,0,1]';
probType = 'tradeOff';
altSolMethod = 'pathwayExclusion';
nSol=5;
[Model_tradeoff_altSols] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType);

