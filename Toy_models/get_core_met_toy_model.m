function model = get_core_met_toy_model()
reactions = {};
for i=1:26
    reactions = [reactions;['r',num2str(i)]];
end
rxnFormula = {'S -> a','b -> c','g -> T2','e -> ii','ii -> T2','f -> jj','jj -> T1','T1 -> T2',...
        'S + m -> n','S -> d','S -> e','S -> f','n -> l','l -> m','l -> jj','k -> jj','T1 -> k','ii -> jj',...
        'd -> h + g','a -> b','h -> T3','c -> T2',' -> S','T3 ->','T2 ->','T1 ->'};
model = createModel(reactions,reactions,rxnFormula);
end