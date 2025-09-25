% this code compares sprintcc stoichiometry based and topology based
% methods for consistency checks
clear
tol=1e-4;

%% creating a toy model with n=1
n=1;
model = topology_toy_model(n);
consistencyType = 'stoichiometry';
ConsReacIDS_stoi_1 = spectraCC(model,tol,consistencyType);

consistencyType = 'topology';
ConsReacIDS_topo_1 = spectraCC(model,tol,consistencyType);

% creating a toy model with n=2
n=2;
model = topology_toy_model(n);
consistencyType = 'stoichiometry';
ConsReacIDS_stoi_2 = spectraCC(model,tol,consistencyType);

consistencyType = 'topology';
ConsReacIDS_topo_2 = spectraCC(model,tol,consistencyType);

% creating a toy model with n=3
n=3;
model = topology_toy_model(n);
consistencyType = 'stoichiometry';
ConsReacIDS_stoi_3 = spectraCC(model,tol,consistencyType);

consistencyType = 'topology';
ConsReacIDS_topo_3 = spectraCC(model,tol,consistencyType);

%% Checking on different toy models
model = get_cc_toy_model_1();
consistencyType = 'stoichiometry';
ConsReacIDS_stoi_11 = spectraCC(model,tol,consistencyType);

consistencyType = 'topology';
ConsReacIDS_topo_11 = spectraCC(model,tol,consistencyType);



model = get_cc_toy_model_2();
consistencyType = 'stoichiometry';
ConsReacIDS_stoi_12 = spectraCC(model,tol,consistencyType);

consistencyType = 'topology';
ConsReacIDS_topo_12 = spectraCC(model,tol,consistencyType);