clear
percents = [10,20,30,40,50];
for i =1:numel(percents)
    GapFillAGORADraftModels(percents(i),'minNetLP',60);
end