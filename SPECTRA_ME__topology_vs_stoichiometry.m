% this code compares sprintgapfiller stoichiometry based and topology based
clear
tol=1e-4;

% creating a toy model with n=1
n=1;
model = topology_toy_model(n);
core = 6; % setting the core reaction to be T1 production

gapFilltype = 'stoichiometry';
[Model_stoichiometry_1,BlockedCoreRxns_stoichiometry_1] = spectraCCME(model,core,tol,gapFilltype);

gapFilltype = 'topology';
[Model_topology_1,BlockedCoreRxns_topology_1] = spectraCCME(model,core,tol,gapFilltype);

% creating a toy model with n=2
n=2;
model = topology_toy_model(n);
core = 6; % setting the core reaction to be T1 production

gapFilltype = 'stoichiometry';
[Model_stoichiometry_2,BlockedCoreRxns_stoichiometry_2] = spectraCCME(model,core,tol,gapFilltype);

gapFilltype = 'topology';
[Model_topology_2,BlockedCoreRxns_topology_2] = spectraCCME(model,core,tol,gapFilltype);

% creating a toy model with n=3
n=3;
model = topology_toy_model(n);
core = 6; % setting the core reaction to be T1 production

gapFilltype = 'stoichiometry';
[Model_stoichiometry_3,BlockedCoreRxns_stoichiometry_3] = spectraCCME(model,core,tol,gapFilltype);

gapFilltype = 'topology';
[Model_topology_3,BlockedCoreRxns_topology_3] = spectraCCME(model,core,tol,gapFilltype);