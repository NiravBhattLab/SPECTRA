clear
load('./consRecon3D')
% initCobraToolbox(false)
changeCobraSolver('gurobi','all')

n = 5000:500:10600; % number of reactions in the inconsistent subnetwork
tol =1e-4; 

nCore = 50:100:5000;
tbl = table('Size',[numel(n)*numel(nCore)*4,5],'VariableTypes',{'string','double','double','double','double'},...
    'VariableNames',{'Algo','InconsModelSize','NoOfCoreRxns','Runtime','LPS'});
k=1;
for i=1:numel(n)
    % reaction IDs of inconsistent subnetwork
    RxnID = sort(randsample(numel(model.rxns),n(i)));
    InconsModel = removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],RxnID)));
    for j=1:numel(nCore)
        core = sort(randsample(numel(InconsModel.rxns),nCore(j)));
        % fastcc+fastcore
        tic
        [a,~,~,LPs1] = fastcc(InconsModel,tol,0);
        m =  removeRxns(InconsModel,InconsModel.rxns(setdiff([1:numel(InconsModel.rxns)],a)));
        cf = find(ismember(m.rxns,InconsModel.rxns(core)));
        [m1,~,LPs2] = fastcore(m,cf,tol);
        time = toc;
        tbl(k,:)={'FC',n(i),nCore(j),time,LPs1+LPs2};
        k=k+1;
        
        % sprintcc+sprintcore
        tic
        [a,LPs1] = sprintcc(InconsModel,tol);
        m =  removeRxns(InconsModel,InconsModel.rxns(setdiff([1:numel(InconsModel.rxns)],a)));
        cf = find(ismember(m.rxns,InconsModel.rxns(core)));
        [m1,LPs2] = sprintcore(m,cf,tol);
        time = toc;
        tbl(k,:)={'SPcc_SPc',n(i),nCore(j),time,LPs1+LPs2};
        k=k+1;
        
        % sprintgapfiller
        tic
        cs = find(ismember([1:numel(InconsModel.rxns)],core));
        [m1,~,LPs] = SprintGapFiller(InconsModel,cs,tol);
        time = toc;
        tbl(k,:)={'SP',n(i),nCore(j),time,LPs};
        k=k+1;
        
        % swiftcc+swiftcore
        tic
        a = swiftcc(InconsModel.S,InconsModel.rev,'gurobi');
        m =  removeRxns(InconsModel,InconsModel.rxns(setdiff([1:numel(InconsModel.rxns)],a)));
        cs = find(ismember(m.rxns,InconsModel.rxns(core)));
        [m1,~,LPs] = swiftcore(m,cs,ones(numel(m.rxns),1),1e-10,0);
        time = toc;
        tbl(k,:)={'SW',n(i),nCore(j),time,LPs+1};
        k=k+1;
    end
end
save('Results_ForSprintGapFillerResults')