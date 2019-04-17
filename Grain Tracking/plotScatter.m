function plotScatter(x, y, label_x, label_y)
% 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3
scatter(x, y,'filled', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3)
xlabel(label_x)
ylabel(label_y)
box on
end