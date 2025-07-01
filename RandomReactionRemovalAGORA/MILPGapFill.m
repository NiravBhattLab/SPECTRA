clear
percents = [10,20];
for i =1:numel(percents)
    GapFillAGORADraftModels(percents(i),'MILP',120);
end