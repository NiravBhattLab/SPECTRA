%% load the cancer-core reactions proportions in each of the models built using fastcore, swiftcore and sprintcore and make the comparisons
load('coverage_fastcore.mat')
load('coverage_swiftcore.mat')
load('coverage_spectra.mat')
spectra_fast_equal=find(coverage_spectra==coverage_fastcore);
spectra_swift_equal=find(coverage_spectra==coverage_swiftcore);
spectra_fast=find(coverage_spectra>coverage_fastcore);
spectra_swift=find(coverage_spectra>coverage_swiftcore);
spectra_fast_less=find(coverage_spectra<coverage_fastcore);
spectra_swift_less=find(coverage_spectra<coverage_swiftcore);
