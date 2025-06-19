clear
load('./consRecon3D')
% initCobraToolbox(false)
changeCobraSolver('gurobi','all')
n = 200:200:10600; % number of core reactions
tol =1e-4; % tolerance level for fastcore and sprintcore

weights = ones(numel(model.rxns),1); %weights for swiftcore

% variables to store the reactions in context-specific model
sprxns={}; sp2rxns={}; sp3rxns={}; swRrxns={}; swWRrxns={}; fcrxns={};
% variabels to store the runtime
sptimes=[]; sp2times=[]; sp3times=[];swRtimes=[]; swWRtimes=[]; fctimes=[];
% variables to store the nLPs
spLPS=[]; sp2LPS=[];swRLPS=[]; swWRLPS=[]; fcLPS=[];
% to store the core reactions obtained by random-sampling at each iteration
cores={};

% This variable stores the information about the models that do not have
% core reactions (or) not consistent
ModelQuality = struct();
ModelQuality.name = {};
ModelQuality.fault = {};
for i=1:numel(n)
    core = sort(randsample(numel(model.rxns),n(i)));
    cores{i} = core;
    core_rxns = model.rxns(core);
    
    % fastcore
    tic
    [m,~,LPs] = fastcore(model,core,tol);
    time = toc;
    fcrxns{i}=m.rxns;
    fctimes(i)=time;
    fcLPS(i)=LPs;
    model_name=['FC_',num2str(i)];
    ModelQuality=checkModelQuality(ModelQuality,m,core_rxns,model_name,tol);
    
    % sprintcore
    tic
    [m,LPs] = sprintcore(model,core,tol);
    time = toc;
    sprxns{i}=m.rxns;
    sptimes(i)=time;
    spLPS(i)=LPs;
    model_name=['SP_',num2str(i)];
    ModelQuality=checkModelQuality(ModelQuality,m,core_rxns,model_name,tol);
    
    % sprintcore2
    tic
    [m,LPs] = sprintcore2(model,core,tol);
    time = toc;
    sp2rxns{i}=m.rxns;
    sp2times(i)=time;
    sp2LPS(i)=LPs;
    model_name=['SP2_',num2str(i)];
    ModelQuality=checkModelQuality(ModelQuality,m,core_rxns,model_name,tol);
    
    % sprintcore3
    tic
    [m,~] = sprintcore(model,core,tol,[],[],[],[],'MILP');
    time = toc;
    sp3rxns{i}=m.rxns;
    sp3times(i)=time;
    model_name=['SP3_',num2str(i)];
    ModelQuality=checkModelQuality(ModelQuality,m,core_rxns,model_name,tol);
    
    % swiftcore with reduction
    tic
    [m,~,LPs] = swiftcore(model,core,weights,1e-10,1);
    time = toc;
    swRrxns{i}=m.rxns;
    swRtimes(i)=time;
    swRLPS(i)=LPs;
    model_name=['SWr_',num2str(i)];
    ModelQuality=checkModelQuality(ModelQuality,m,core_rxns,model_name,tol);
    
    % swiftcore without reduction
    tic
    [m,~,LPs] = swiftcore(model,core,weights,1e-10,0);
    time = toc;
    swWRrxns{i}=m.rxns;
    swWRtimes(i)=time;
    swWRLPS(i)=LPs;
    model_name=['SWwr_',num2str(i)];
    ModelQuality=checkModelQuality(ModelQuality,m,core_rxns,model_name,tol);
end
save('Results_SprintCoreRuntime')