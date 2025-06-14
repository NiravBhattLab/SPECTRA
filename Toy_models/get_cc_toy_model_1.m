function model = get_cc_toy_model_1()
    reactions = {'r1','r2','r3','r4','r5','r6','r7','r8'};
    rxnFormula = {'-> A','A -> B','B -> C','B -> D','D ->','A -> E','E -> F','F -> D'};
    model = createModel(reactions,reactions,rxnFormula);
end