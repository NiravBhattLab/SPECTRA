clear
load('consRecon3D') % loading the consistent Recon3D model
changeCobraSolver('gurobi','all')
model.rev = model.lb<0;
n=20:20:10600;
n_lp7=[];
n_lpf=[];
tol = 1e-4;
for i =1:numel(n)
    rxn_ids = sort(randsample(numel(model.rxns),n(i)));
    m = removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],rxn_ids)));
   [flux,~] = forwardcc(m,true(numel(m.rxns),1),tol);
   n_lpf(i) = sum(abs(flux)>=tol);
   [flux,~] = fc_lp7(m,true(numel(m.rxns),1),tol);
   n_lp7(i) = sum(abs(flux)>=tol);
end
figure()
plot(n,n_lp7,'.r','MarkerSize',5)
hold on
plot(n,n_lpf,'.g','MarkerSize',5)
xlabel('Number of reactions in the model','fontweight','bold','fontsize',20)
ylabel({'Number of reactions identified to be'; 'consistent in first iteration'},'fontweight','bold','fontsize',20)
getLegend(2,{'\bf LP7','\bf LPforwardCC'},{'r';'g'},15,'northwest')
set(gca,'fontweight','bold')

function getLegend(n,leg,colors,fontsize,lcn)

% Now Create bar chart with nan so it won't show:
b = bar(nan(n,n));
set(b,{'FaceColor'},colors)
lgd = legend(b,leg,'Location',lcn);
lgd.FontSize=fontsize;
end