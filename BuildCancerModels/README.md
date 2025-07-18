# Folder description

**Requirements to reproduce the codes**
1) The cancer cell line transcriptomics data of 1479 cell lines from the Cancer Cell Line Encyclopedia (CCLE) (Refer: https://depmap.org/portal/)


1. DepMapModels_Recon (To build the required context-specific models using distict MeMs)
- PrepareDepMapData.m to extract the transcriptomics data from the CCLE database into .mat structures. The data is stored in DepMapData.mat
- modGenesIDConv.m to convert the gene ids provided in the UpdatedRecon3D.mat in the form of Recon 3D id format into official gene names. The updated model is consRecon3DGeneSymbol.mat
- BuildDepMapModels.m to build the models using fastcore, swiftcore and sprintcore using consRecon3DGeneSymbol.mat and DepMapData.mat
- GiniReactionImportance.m to get the list of core reactions from the transcriptomics data

2. coverage_core
- Cancerrxnscoverage.m to calculate the proportion of generic core cancer reactions in each of the models built using fastcore, swiftcore and sprintcore
- core_coverage_find.m to compare the results core cancer reactions coverage results of sprintcore derived models with the models built using fastcore and swiftcore

3. fva
- Excel files that has the details of the min and max flux obtained using FVA for the 20 models in the study.

4. flux_results
- models_name.m contains the model names for cancerous and the corresponding non-cancerous contexts
- readfva.m extracts the flux variability results from the excel for all models and stores in .mat format
- fsranalysis_ccle_sprint.m to calculate the flux span ratios and flux enrichment ratios across the 20 models
