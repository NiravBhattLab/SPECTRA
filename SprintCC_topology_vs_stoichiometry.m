% this code compares sprintcc stoichiometry based and topology based
% methods for consistency checks
clear
tol=1e-4;

%% creating a toy model with n=1
n=1;
model = topology_toy_model(n);
consistencyType = 'stoichiometry';
ConsReacIDS_stoi_1 = sprintcc(model,tol,consistencyType);

consistencyType = 'topology';
ConsReacIDS_topo_1 = sprintcc(model,tol,consistencyType);

% creating a toy model with n=2
n=2;
model = topology_toy_model(n);
consistencyType = 'stoichiometry';
ConsReacIDS_stoi_2 = sprintcc(model,tol,consistencyType);

consistencyType = 'topology';
ConsReacIDS_topo_2 = sprintcc(model,tol,consistencyType);

% creating a toy model with n=3
n=3;
model = topology_toy_model(n);
consistencyType = 'stoichiometry';
ConsReacIDS_stoi_3 = sprintcc(model,tol,consistencyType);

consistencyType = 'topology';
ConsReacIDS_topo_3 = sprintcc(model,tol,consistencyType);

%% Checking on a different toy model
model = get_cc_toy_model_1();
consistencyType = 'stoichiometry';
ConsReacIDS_stoi_11 = sprintcc(model,tol,consistencyType);

consistencyType = 'topology';
ConsReacIDS_topo_11 = sprintcc(model,tol,consistencyType);

model = get_cc_toy_model_2();
consistencyType = 'stoichiometry';
ConsReacIDS_stoi_12 = sprintcc(model,tol,consistencyType);

consistencyType = 'topology';
ConsReacIDS_topo_12 = sprintcc(model,tol,consistencyType);