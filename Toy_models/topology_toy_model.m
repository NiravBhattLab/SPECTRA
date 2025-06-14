function model = topology_toy_model(n)
% this code generates a toy model as shown in Figure 2a) of the article
% "Meneco, a Topology-Based Gap-Filling Tool Applicable to Degraded
%  Genome-Wide Metabolic Networks"

reactions = {'r1','r2','r3','r4','r5','r6'};
rxnFormula = {'-> S', ['S -> b + ',num2str(n),' a'],...
    'b + a -> d','a -> c','d + c -> T1','T1 ->'};

model = createModel(reactions,reactions,rxnFormula);
end