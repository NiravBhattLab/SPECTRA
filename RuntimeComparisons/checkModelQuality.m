function ModelQuality=checkModelQuality(ModelQuality,m,core_rxns,model_name,tol)
% This function checks if the context-specific model has all the core
% reactions in it and check if the model is consistent
rev = m.lb<0;
a = swiftcc(m.S,rev,'gurobi');
if numel(a)~=numel(m.rxns)
    ModelQuality.name = [ModelQuality.name;model_name];
    ModelQuality.fault = [ModelQuality.fault;'inconsistent'];
end
if sum(ismember(core_rxns,m.rxns))~=numel(core_rxns)
    ModelQuality.name = [ModelQuality.name;model_name];
    ModelQuality.fault = [ModelQuality.fault;'Missing core rxns'];
end
end