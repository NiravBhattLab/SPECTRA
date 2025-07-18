clear
percents = [10];
for i =1:numel(percents)
    GapFillAGORADraftModels(percents(i),'MILP',180);
end