rxns = {'r1','r2','r3','r4','r5','r6','r7','r8','r9'};
rxnform = {'-> A','A -> B','B ->','-> A','B ->','-> A','B ->','-> A','B ->'};
model = createModel(rxns,rxns,rxnform);
[Model2,BlockedCoreRxns,flux,LPs] = SprintGapFiller(model,2,1e-4,[],20,'pathwayExclusion','MILP');
[Model1,LPS] = sprintcore(model,2,1e-4,[],20,'pathwayExclusion','MILP');


temp1 = [];
for i=1:numel(Model1)
    temp = Model1{i};
    temp1 = [temp1;strjoin(sort(string(temp.rxns)))];
end
temp2 = [];
for i=1:numel(Model2)
    temp = Model2{i};
    temp2 = [temp2;strjoin(sort(string(temp.rxns)))];
end