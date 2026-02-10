clear
model = get_rand_wts_toy_model;
core = [2,4,5];
csm_fastcore = fastcore(model,core,1e-4);
csm_spectra = spectraME(model,core,1e-4);