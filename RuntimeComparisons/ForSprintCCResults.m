clear
load('consRecon3D')
% initCobraToolbox(false)
changeCobraSolver('gurobi','all')
model.rev = model.lb<0;
tol=1e-4;
n=20:20:10600;
t_fcc=[];t_spcc=[];t_swcc=[];
a_fcc={};a_spcc={};a_swcc={};
LPS_fcc=[];LPS_spcc=[];
cores={};
for i=1:numel(n)
    rxn_ids = sort(randsample(numel(model.rxns),n(i)));
    cores{i}=rxn_ids;
    m = removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],rxn_ids)));
    
    tic
    [a1,~,~,LPs] = fastcc(m,tol);
    t_fcc(i)=toc;
    a_fcc{i} = a1;
    LPS_fcc(i) = LPs;
    
    tic
    a2 = swiftcc(m.S,m.rev,'gurobi');
    t_swcc(i)=toc;
    a_swcc{i} = a2;
    
    tic
    [a3,LPs] = sprintcc(m,tol);
    t_spcc(i)=toc;
    a_spcc{i} = a3;
    LPS_spcc(i)=LPs;
end
tbl = table('Size',[530*3,4],'VariableTypes',{'string','double','double','double'},...
    'VariableNames',{'Algo','SizeOfModel','Runtime','LPS'});
col1 = [repmat({'FastCC'},530,1);repmat({'SwiftCC'},530,1);repmat({'SprintCC'},530,1)];
col2 = repmat(n',3,1);
col3 = [t_fcc';t_swcc';t_spcc'];
col4 = [LPS_fcc';ones(530,1);LPS_spcc'];
tbl.Algo =col1;
tbl.SizeOfModel=col2;
tbl.Runtime = col3;
tbl.LPS = col4;
