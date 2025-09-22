function model = three_pathway_toy_model()
% this code generates a toy model with three alternate pathways to 
% produce the biomass metabolite"

reactions = arrayfun(@(x)['r',num2str(x)],1:11,'UniformOutput',false);
rxnFormula = {' -> a','a -> b','b -> c','c -> d','d -> ','-> e','e -> f','f -> 0.5 d',' -> g','g -> h','h -> 0.5 d'};
model = createModel(reactions,reactions,rxnFormula);
end