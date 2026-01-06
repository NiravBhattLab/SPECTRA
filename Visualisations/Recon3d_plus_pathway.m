clear
load('..\Recon3D+\UpdatedRecon3D.mat')
r3d_plus = model;
load('..\LP7_vs_LPf\consRecon3D.mat')
r3d = model;

new_rxns = setdiff(r3d_plus.rxns,r3d.rxns);
ss = r3d_plus.subSystems(ismember(r3d_plus.rxns,new_rxns));
ss = [ss{:}]';

[uniqNames, ~, idx] = unique(ss);
counts = accumarray(idx, 1);

[countsSorted, order] = sort(counts, 'descend');
uniqNamesSorted = uniqNames(order);

N = min(20, numel(countsSorted));
countsTop = flip(countsSorted(1:N));
namesTop  = flip(uniqNamesSorted(1:N));
namesTop{N} = {'Transport'; 'endoplasmic reticular'};
namesTop{N-1} = {'Exchange/demand'; 'reaction'};
namesTop{N-2} = {'Glycerophospholipid'; 'metabolism'};
maxLen = 10;
lineGap = 0.5;
FontSize = 15;
figure
b = barh(countsTop, 'FaceColor', '#008080');
set(gca, 'YTickLabel', {}, 'FontSize', FontSize);

for i = 1:N
    if ismember(i,[N-2:N])
        text(countsTop(i) + 1, i + lineGap/2, namesTop{i}{1}, ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','middle', ...
            'FontSize',FontSize);
        text(countsTop(i) + 1, i - lineGap/2, namesTop{i}{2}, ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','middle', ...
            'FontSize',FontSize);
    else
        text(countsTop(i) + 1, i, ...
             namesTop{i}, ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'middle', ...
             'FontSize', FontSize);
    end
end
xlabel('Number of reactions','FontWeight','bold','FontSize',20)
ylabel('Subsystems','FontWeight','bold','FontSize',20)
xticks([0:30:150])
xticklabels([0:30:150])
title('Top 20 subsystems of newly added reactions in Recon3D+')
xlim([0,150])
ax = gca;
ax.Box = 'off';

set(gcf, 'Renderer', 'painters');  % Vector graphics
print('Top_20_recon3d_pathways.png', '-dpng', '-r300');