contexts={'Brain'
'Kidney'
'Liver'
'Lung'
'Eye'
'Pancreas'
'Prostate'
'Breast'
'Ovary'
'Biliary Tract'};
%%
p=dir('../fva');
p = {p(1:end).name}';
p=p(find(contains(p,'ACH-')));

run models_name.m;
load('cells_recon_2025.mat');
for j =1:length(contexts)
    model_normal_new=models_normal{j,1};
    model_diseased_new=models_cancer{j,1};
    %%
    fsr_files = split(p,'.xlsx');
    fsr_files = append(fsr_files,'.mat');
    %%
    matches(fsr_files(:,1),model_normal_new);
    find(ans==1);
    cellsinclude_1=cells_all{1,ans};

    %%
    matches(fsr_files(:,1),model_diseased_new);
    find(ans==1);
    cellsinclude_2=cells_all{1,ans};
    %%
    commonrxns_all=intersect(cellsinclude_1(:,1),cellsinclude_2(:,1));
    %%
    rxns_normal=cellsinclude_1(:,1);
    fsrcase_normal=cell2mat(cellsinclude_1(:,2:3));
    %%
    rxns_cancer=cellsinclude_2(:,1);
    fsrcase_cancer=cell2mat(cellsinclude_2(:,2:3));
    %%
    minFlux0=fsrcase_normal(:,1);
    minFlux=fsrcase_cancer(:,1);
    maxFlux0=fsrcase_normal(:,2);
    maxFlux=fsrcase_cancer(:,2);
    %%
    minFlux_normal=(minFlux0);
    maxFlux_normal=(maxFlux0);
    minFlux_diseased=(minFlux);
    maxFlux_diseased=(maxFlux);
    model_normal=load(['../DepMapModels_Recon/sprintcore/',model_normal_new]);
    model_diseased=load(['../DepMapModels_Recon/sprintcore/',model_diseased_new]);
   %%
    exchanges_all=commonrxns_all(contains(commonrxns_all,'EX_'));
    sinks_all=commonrxns_all(contains(commonrxns_all,'sink_'));
    demands_all=commonrxns_all(contains(commonrxns_all,'DM_'));
    %%
    commonrxns = setdiff(commonrxns_all,[exchanges_all;sinks_all;demands_all]);
    rxnids_normal=findRxnIDs(model_normal.m,commonrxns);
    rxnids_diseased=findRxnIDs(model_diseased.m,commonrxns);
    maxFlux_normal_met=maxFlux_normal(rxnids_normal);
    minFlux_normal_met=minFlux_normal(rxnids_normal);
    maxFlux_diseased_met=maxFlux_diseased(rxnids_diseased);
    minFlux_diseased_met=minFlux_diseased(rxnids_diseased);
    %%
    fluxspanrationew0=[];
    for i = 1:length(commonrxns)
	    %if maxFlux_diseased_met(i) > minFlux_diseased_met(i) && maxFlux_normal_met(i) > minFlux_normal_met(i); 
         if maxFlux_diseased_met(i) ~= maxFlux_normal_met(i) && minFlux_diseased_met(i) <= minFlux_normal_met(i);
	        fluxspanrationew0(i)=(maxFlux_diseased_met(i) - minFlux_diseased_met(i))./(maxFlux_normal_met(i) - minFlux_normal_met(i));
         elseif maxFlux_diseased_met(i) == maxFlux_normal_met(i) && minFlux_diseased_met(i) == minFlux_normal_met(i);
            fluxspanrationew0(i)=0;
	    %else fluxspanrationew0(i)=(maxFlux_diseased(i) - minFlux_diseased(i))./maxFlux_normal(i);
        end
    end
    fluxspanrationew=[fluxspanrationew0]';
    %% reactions
    Downregulated_disease =commonrxns(find((fluxspanrationew~=0 & isfinite(fluxspanrationew) & ~isnan(fluxspanrationew) & fluxspanrationew<=0.5 & fluxspanrationew>=0.01)));
    Upregulated_disease = commonrxns(fluxspanrationew~=0 & isfinite(fluxspanrationew) & ~isnan(fluxspanrationew) & fluxspanrationew>=2 & fluxspanrationew<20000.0);
    %% fsr values
    fsr_Downregulated_disease=fluxspanrationew(find((fluxspanrationew~=0 & isfinite(fluxspanrationew) & ~isnan(fluxspanrationew) & fluxspanrationew<=0.5 & fluxspanrationew>=0.01)));
    fsr_Upregulated_disease=fluxspanrationew(find((fluxspanrationew~=0 & isfinite(fluxspanrationew) & ~isnan(fluxspanrationew) & fluxspanrationew>=2& fluxspanrationew<20000.0)));
    %% upregulated reactions
    Upregulated_disease_ids = findRxnIDs(model_normal.m,Upregulated_disease);
    Upregulated_disease_names=model_normal.m.rxnNames(Upregulated_disease_ids);
    Upregulated_disease_subsystems=model_normal.m.subSystems(Upregulated_disease_ids);
    Upregulated_disease_formulae=printRxnFormula(model_normal.m,Upregulated_disease);
    Upregulated_disease_genes=model_normal.m.grRules(Upregulated_disease_ids);
    display(length(Upregulated_disease_subsystems))
    Upregulated_disease_details=table(Upregulated_disease,fsr_Upregulated_disease,Upregulated_disease_names,Upregulated_disease_formulae,Upregulated_disease_genes,Upregulated_disease_subsystems);
    %% downregulated reactions
    Down_disease_ids = findRxnIDs(model_normal.m,Downregulated_disease);
    Down_disease_names=model_normal.m.rxnNames(Down_disease_ids);
    Down_disease_subsystems=model_normal.m.subSystems(Down_disease_ids);
    Down_disease_formulae=printRxnFormula(model_normal.m,Downregulated_disease);
    Down_disease_genes=model_normal.m.grRules(Down_disease_ids);
    Down_disease_details=table(Downregulated_disease,fsr_Downregulated_disease,Down_disease_names,Down_disease_formulae,Down_disease_genes,Down_disease_subsystems);
    %%
    regulated_diseases_fsr=[table2cell(Upregulated_disease_details(:,2));table2cell(Down_disease_details(:,2))];
    regulated_diseases=[table2cell(Upregulated_disease_details(:,1));table2cell(Down_disease_details(:,1))];
    regulated_diseases_names=[table2cell(Upregulated_disease_details(:,3));table2cell(Down_disease_details(:,3))];
    regulated_diseases_subsystems=[table2cell(Upregulated_disease_details(:,4));table2cell(Down_disease_details(:,4))];
    table2cell(Down_disease_details)
    downrxns=ans;
    table2cell(Upregulated_disease_details);
    uprxns=ans;
    writetable(Upregulated_disease_details,'fsr_up_2025_checked_2.xlsx','Sheet',contexts{j,1});
    writetable(Down_disease_details,'fsr_down_2025_checked_2.xlsx','Sheet',contexts{j,1});
    clear fluxspanrationew0 fluxspanrationew
    diffids=findRxnIDs(model_diseased.m,[Downregulated_disease;Upregulated_disease]);
    resultCell = FEA(model_diseased.m,diffids,'subSystems');
    featable=table(resultCell)
    writetable(featable,'fea_2025_checked_2.xlsx','Sheet',contexts{j,1});
end
