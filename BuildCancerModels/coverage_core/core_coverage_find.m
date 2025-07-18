%% load the cancer-core reactions proportions in each of the models built using fastcore, swiftcore and sprintcore and make the comparisons
load('coverage_fastcore.mat')
load('coverage_swiftcore.mat')
load('coverage_sprintcore.mat')
sprint_fast_equal=find(coverage_sprintcore==coverage_fastcore);
sprint_swift_equal=find(coverage_sprintcore==coverage_swiftcore);
sprint_fast=find(coverage_sprintcore>coverage_fastcore);
sprint_swift=find(coverage_sprintcore>coverage_swiftcore);
sprint_fast_less=find(coverage_sprintcore<coverage_fastcore);
sprint_swift_less=find(coverage_sprintcore<coverage_swiftcore);
