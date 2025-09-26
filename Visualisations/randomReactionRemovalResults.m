clear
load('./Results/getResultsOnGapfill30_LP_MILP_DC_2.mat')
acc_LP = (tp_LP+tn_LP)./(tp_LP+tn_LP+fp_LP+fn_LP);
acc_MILP = (tp_MILP+tn_MILP)./(tp_MILP+tn_MILP+fp_MILP+fn_MILP);
acc_DC = (tp_DC+tn_DC)./(tp_DC+tn_DC+fp_DC+fn_DC);

figure()
titles ={'LP';'MILP';'DC'};
boxchart([acc_LP',acc_MILP',acc_DC'])
set(gca,'XTickLabel',titles,'fontsize',12,'FontWeight','bold')
xlabel('Algorithm','fontweight','bold','fontsize',20)
ylabel('Accuracy','fontweight','bold','fontsize',20)

figure()
titles ={'LP';'MILP';'DC'};
boxchart([f1_LP',f1_MILP',f1_DC'])
set(gca,'XTickLabel',titles,'fontsize',12,'FontWeight','bold')
xlabel('Algorithm','fontweight','bold','fontsize',20)
ylabel('f1 score','fontweight','bold','fontsize',20)

figure()
boxchart([tp_LP',tp_MILP',tp_DC'])
set(gca,'XTickLabel',titles,'fontsize',12,'FontWeight','bold')
xlabel('Algorithm','fontweight','bold','fontsize',20)
ylabel('TP','fontweight','bold','fontsize',20)

figure()
boxchart([tn_LP',tn_MILP',tn_DC'])
set(gca,'XTickLabel',titles,'fontsize',12,'FontWeight','bold')
xlabel('Algorithm','fontweight','bold','fontsize',20)
ylabel('TN','fontweight','bold','fontsize',20)

figure()
boxchart([fn_LP',fn_MILP',fn_DC'])
set(gca,'XTickLabel',titles,'fontsize',12,'FontWeight','bold')
xlabel('Algorithm','fontweight','bold','fontsize',20)
ylabel('FN','fontweight','bold','fontsize',20)

figure()
boxchart([fp_LP',fp_MILP',fp_DC'])
set(gca,'XTickLabel',titles,'fontsize',12,'FontWeight','bold')
xlabel('Algorithm','fontweight','bold','fontsize',20)
ylabel('FP','fontweight','bold','fontsize',20)

figure()
boxchart([pre_LP',pre_MILP',pre_DC'])
set(gca,'XTickLabel',titles,'fontsize',12,'FontWeight','bold')
xlabel('Algorithm','fontweight','bold','fontsize',20)
ylabel('Precision','fontweight','bold','fontsize',20)

figure()
boxchart([rec_LP',rec_MILP',rec_DC'])
set(gca,'XTickLabel',titles,'fontsize',12,'FontWeight','bold')
xlabel('Algorithm','fontweight','bold','fontsize',20)
ylabel('Recall','fontweight','bold','fontsize',20)