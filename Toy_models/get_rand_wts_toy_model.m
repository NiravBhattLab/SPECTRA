function model = get_rand_wts_toy_model()
  reactions = {'r1','r2','r3','r4','r5','r6','r7','r8','r9'};
  rxnFormula = {'-> A','A <=> B','B -> C','C -> D','E <=> D','E ->','B ->','D <=>','-> C'};
  model = createModel(reactions,reactions,rxnFormula);
end