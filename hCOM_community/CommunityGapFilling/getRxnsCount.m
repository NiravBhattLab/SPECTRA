% to get the count of reactions in the gap filled models
clear
p = dir('./CarvemeCommunityGapfillingLP');
p = {p(3:end).name}';

single_lp=[];
single_milp=[];
comm_lp=[];
comm_milp=[];
for i = 1:numel(p)
    load(['./CarvemeCommunityGapfillingLP/',p{i}])
    comm_lp=[comm_lp;numel(Model.rxns)];

    load(['./CarvemeCommunityGapfillingMILP/',p{i}])
    comm_milp = [comm_milp;numel(Model.rxns)];

    load(['./CarvemeIndividualGapfillingLP/',p{i}])
    single_lp = [single_lp; numel(Model.rxns)];

    load(['./CarvemeIndividualGapfillingMILP/',p{i}])
    single_milp = [single_milp; numel(Model.rxns)];
end
save('Results_getRxnsCount')