clear
load('./Results/getResultsOnGapfill30_LP_MILP_DC_2.mat')
acc_LP = (tp_LP+tn_LP)./(tp_LP+tn_LP+fp_LP+fn_LP);
acc_MILP = (tp_MILP+tn_MILP)./(tp_MILP+tn_MILP+fp_MILP+fn_MILP);
acc_DC = (tp_DC+tn_DC)./(tp_DC+tn_DC+fp_DC+fn_DC);

titles = {'minNetLP','minNetMILP','minNetDC'};
% helper inline function to convert hex -> RGB
hex2rgb = @(hex) sscanf(hex(2:end),'%2x%2x%2x',[1 3])/255;
% define colors
c1 = hex2rgb('#008080'); 
c2 = hex2rgb('#004080'); 
c3 = hex2rgb('#80c080');


sigLabel = @(p) ...
    ternary(p < 0.05, '*', 'n.s.');



%%
figure()
subplot(2,2,1)
boxchart(ones(size(tp_LP')),    tp_LP',   'BoxFaceColor', c1, 'MarkerColor', c1, 'LineWidth', 2);
hold on
boxchart(2*ones(size(tp_MILP')), tp_MILP','BoxFaceColor', c2, 'MarkerColor', c2, 'LineWidth', 2);
boxchart(3*ones(size(tp_DC')),   tp_DC',  'BoxFaceColor', c3, 'MarkerColor', c3, 'LineWidth', 2);
hold off
set(gca,'XTick',1:3,'XTickLabel',titles,'FontSize',20,'FontWeight','bold')
xlabel('Algorithm','fontweight','bold','fontsize',20)
ylabel('True Positives','fontweight','bold','fontsize',20)

hold on
[p_MILP_LP, ~] = signrank(tp_MILP, tp_LP, 'tail','right');
[p_MILP_DC,~] = signrank(tp_MILP, tp_DC, 'tail','right');
% x positions used in your plot
x1 = 1;
x2 = 2;
x3 = 3;
% height for bars
y = max([tp_LP, tp_MILP, tp_DC]) + 0.05;
h = 0.03;
% LP vs MILP
plot([x1 x1 x2 x2], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x1 x2]), y+2*h, sigLabel(p_MILP_LP), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
% MILP vs DC
y = y + 0.08;
plot([x2 x2 x3 x3], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x2 x3]), y+2*h, sigLabel(p_MILP_DC), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
hold off

%%
subplot(2,2,2)
boxchart(ones(size(tn_LP')),    tn_LP',   'BoxFaceColor', c1, 'MarkerColor', c1, 'LineWidth', 2);
hold on
boxchart(2*ones(size(tn_MILP')), tn_MILP','BoxFaceColor', c2, 'MarkerColor', c2, 'LineWidth', 2);
boxchart(3*ones(size(tn_DC')),   tn_DC',  'BoxFaceColor', c3, 'MarkerColor', c3, 'LineWidth', 2);
hold off
set(gca,'XTick',1:3,'XTickLabel',titles,'FontSize',20,'FontWeight','bold')
xlabel('Algorithm','fontweight','bold','fontsize',20)
ylabel('True Negatives','fontweight','bold','fontsize',20)

hold on
[p_MILP_LP, ~] = signrank(tn_MILP, tn_LP, 'tail','right');
[p_MILP_DC,~] = signrank(tn_MILP, tn_DC, 'tail','right');
% x positions used in your plot 
x1 = 1;
x2 = 2;
x3 = 3;
% height for bars
y = max([tn_LP, tn_MILP, tn_DC]) + 0.05;
h = 0.03;
% LP vs MILP
plot([x1 x1 x2 x2], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x1 x2]), y+2*h, sigLabel(p_MILP_LP), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
% MILP vs DC
y = y + 0.08;
plot([x2 x2 x3 x3], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x2 x3]), y+2*h, sigLabel(p_MILP_DC), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
hold off

%%
subplot(2,2,3)
boxchart(ones(size(fn_LP')),    fn_LP',   'BoxFaceColor', c1, 'MarkerColor', c1, 'LineWidth', 2);
hold on
boxchart(2*ones(size(fn_MILP')), fn_MILP','BoxFaceColor', c2, 'MarkerColor', c2, 'LineWidth', 2);
boxchart(3*ones(size(fn_DC')),   fn_DC',  'BoxFaceColor', c3, 'MarkerColor', c3, 'LineWidth', 2);
hold off
set(gca,'XTick',1:3,'XTickLabel',titles,'FontSize',20,'FontWeight','bold')
xlabel('Algorithm','fontweight','bold','fontsize',20)
ylabel('False Negatives','fontweight','bold','fontsize',20)

hold on
[p_MILP_LP, ~] = signrank(fn_MILP, fn_LP, 'tail','left');
[p_MILP_DC,~] = signrank(fn_MILP, fn_DC, 'tail','left');
% x positions used in your plot
x1 = 1;
x2 = 2;
x3 = 3;
% height for bars
y = max([fn_LP, fn_MILP, fn_DC]) + 0.05;
h = 0.03;
% LP vs MILP
plot([x1 x1 x2 x2], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x1 x2]), y+2*h, sigLabel(p_MILP_LP), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
% MILP vs DC
y = y + 0.08;
plot([x2 x2 x3 x3], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x2 x3]), y+2*h, sigLabel(p_MILP_DC), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
hold off

%%
subplot(2,2,4)
boxchart(ones(size(fp_LP')),    fp_LP',   'BoxFaceColor', c1, 'MarkerColor', c1, 'LineWidth', 2);
hold on
boxchart(2*ones(size(fp_MILP')), fp_MILP','BoxFaceColor', c2, 'MarkerColor', c2, 'LineWidth', 2);
boxchart(3*ones(size(fp_DC')),   fp_DC',  'BoxFaceColor', c3, 'MarkerColor', c3, 'LineWidth', 2);
hold off
set(gca,'XTick',1:3,'XTickLabel',titles,'FontSize',20,'FontWeight','bold')
xlabel('Algorithm','fontweight','bold','fontsize',20)
ylabel('False Positives','fontweight','bold','fontsize',20)

hold on
[p_MILP_LP, ~] = signrank(fp_MILP, fp_LP, 'tail','left');
[p_MILP_DC,~] = signrank(fp_MILP, fp_DC, 'tail','left');
% x positions used in your plot
x1 = 1;
x2 = 2;
x3 = 3;
% height for bars
y = max([fp_LP, fp_MILP, fp_DC]) + 0.05;
h = 0.03;
% LP vs MILP
plot([x1 x1 x2 x2], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x1 x2]), y+2*h, sigLabel(p_MILP_LP), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
% MILP vs DC
y = y + 0.08;
plot([x2 x2 x3 x3], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x2 x3]), y+2*h, sigLabel(p_MILP_DC), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
hold off

%%
figure()
% Use tiledlayout to control spacing
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')  % reduces space between subplots

%% First subplot: Accuracy
ax1 = nexttile;
boxchart(0.4*ones(size(acc_LP')),    acc_LP',   'BoxFaceColor', c1, 'MarkerColor', c1, 'LineWidth', 2);
hold on
boxchart(ones(size(acc_MILP')), acc_MILP','BoxFaceColor', c2, 'MarkerColor', c2, 'LineWidth', 2);
boxchart(1.6*ones(size(acc_DC')),   acc_DC',  'BoxFaceColor', c3, 'MarkerColor', c3, 'LineWidth', 2);
hold off
set(gca,'LineWidth',2) 
set(ax1,'YTick',[0,0.5,1],'XTick',[],'FontSize',20,'FontWeight','bold')
xlim([0.1, 1.9])
ylim([0, 1.2])
ylabel('Accuracy','FontWeight','bold','FontSize',30)

hold on
[p_MILP_LP, ~] = signrank(acc_MILP, acc_LP, 'tail','right');
[p_MILP_DC,~] = signrank(acc_MILP, acc_DC, 'tail','right');
% x positions used in your plot
x1 = 0.4;
x2 = 1.0;
x3 = 1.6;
% height for bars
y = 1.05;
h = 0.03;
% LP vs MILP
plot([x1 x1 x2 x2], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x1 x2]), y+2*h, sigLabel(p_MILP_LP), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
% MILP vs DC
y = y + 0.08;
plot([x2 x2 x3 x3], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x2 x3]), y+2*h, sigLabel(p_MILP_DC), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
hold off


%% Second subplot: F1 score
ax2 = nexttile;
boxchart(0.4*ones(size(f1_LP')),    f1_LP',   'BoxFaceColor', c1, 'MarkerColor', c1, 'LineWidth', 2);
hold on
boxchart(ones(size(f1_MILP')), f1_MILP','BoxFaceColor', c2, 'MarkerColor', c2, 'LineWidth', 2);
boxchart(1.6*ones(size(f1_DC')),   f1_DC',  'BoxFaceColor', c3, 'MarkerColor', c3, 'LineWidth', 2);
hold off
set(gca,'LineWidth',2) 
set(ax2,'YTick',[0,0.5,1],'XTick',[],'FontSize',15,'FontWeight','bold')
xlim([0.1, 1.9])
ylim([0, 1.2])
ylabel('F1 score','FontWeight','bold','FontSize',30)


hold on
[p_MILP_LP, ~] = signrank(f1_MILP, f1_LP, 'tail','right');
[p_MILP_DC,~] = signrank(f1_MILP, f1_DC, 'tail','right');
% x positions used in your plot
x1 = 0.4;
x2 = 1.0;
x3 = 1.6;
% height for bars
y = 1.05;
h = 0.03;
% LP vs MILP
plot([x1 x1 x2 x2], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x1 x2]), y+2*h, sigLabel(p_MILP_LP), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
% MILP vs DC
y = y + 0.08;
plot([x2 x2 x3 x3], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x2 x3]), y+2*h, sigLabel(p_MILP_DC), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
hold off


%%
figure()
% Use tiledlayout to control spacing
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')  % reduces space between subplots

%% First subplot
ax1 = nexttile;
boxchart(0.4*ones(size(pre_LP')),    pre_LP',   'BoxFaceColor', c1, 'MarkerColor', c1, 'LineWidth', 2);
hold on
boxchart(ones(size(pre_MILP')), pre_MILP','BoxFaceColor', c2, 'MarkerColor', c2, 'LineWidth', 2);
boxchart(1.6*ones(size(pre_DC')),   pre_DC',  'BoxFaceColor', c3, 'MarkerColor', c3, 'LineWidth', 2);
hold off
set(gca,'LineWidth',2) 
set(ax1,'YTick',[0,0.5,1],'XTick',[],'FontSize',15,'FontWeight','bold')
xlim([0.1, 1.9])
ylim([0, 1.2])
ylabel('Precision','fontweight','bold','fontsize',30)

hold on
[p_MILP_LP, ~] = signrank(pre_MILP, pre_LP, 'tail','right');
[p_MILP_DC,~] = signrank(pre_MILP, pre_DC, 'tail','right');
% x positions used in your plot
x1 = 0.4;
x2 = 1.0;
x3 = 1.6;
% height for bars
y = 1.05;
h = 0.03;
% LP vs MILP
plot([x1 x1 x2 x2], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x1 x2]), y+2*h, sigLabel(p_MILP_LP), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
% MILP vs DC
y = y + 0.08;
plot([x2 x2 x3 x3], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x2 x3]), y+2*h, sigLabel(p_MILP_DC), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
hold off


%% Second subplot
ax2 = nexttile;
boxchart(0.4*ones(size(rec_LP')),    rec_LP',   'BoxFaceColor', c1, 'MarkerColor', c1, 'LineWidth', 2);
hold on
boxchart(ones(size(rec_MILP')), rec_MILP','BoxFaceColor', c2, 'MarkerColor', c2, 'LineWidth', 2);
boxchart(1.6*ones(size(rec_DC')),   rec_DC',  'BoxFaceColor', c3, 'MarkerColor', c3, 'LineWidth', 2);
hold off
set(gca,'LineWidth',2) 
set(ax2,'YTick',[0,0.5,1],'XTick',[0.4,1,1.6],'XTickLabel',titles,'FontSize',16,'FontWeight','bold')
xlim([0.1, 1.9])
ylim([0, 1.2])
xlabel('Algorithm','fontweight','bold','fontsize',30)
ylabel('Recall','fontweight','bold','fontsize',30)

hold on
[p_MILP_LP, ~] = signrank(rec_MILP, rec_LP, 'tail','right');
[p_MILP_DC,~] = signrank(rec_MILP, rec_DC, 'tail','right');
% x positions used in your plot
x1 = 0.4;
x2 = 1.0;
x3 = 1.6;
% height for bars
y = 1.05;
h = 0.03;
% LP vs MILP
plot([x1 x1 x2 x2], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x1 x2]), y+2*h, sigLabel(p_MILP_LP), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
% MILP vs DC
y = y + 0.08;
plot([x2 x2 x3 x3], [y y+h y+h y], 'k', 'LineWidth', 2)
text(mean([x2 x3]), y+2*h, sigLabel(p_MILP_DC), ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold')
hold off

%%
function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end