
set(0,'defaultAxesLabelFontSize',1.2)
set(0,'defaultAxesFontSize',17)

ylim_bin = [-1e4, 1e4];
figure(1)
plotData_1 = [FItgCurvs_An4(faceCorresp(:,1),2), FA_diff];
mask = (plotData_1(:,2) > ylim_bin(1) & plotData_1(:,2) < ylim_bin(2));
plotData_1 = plotData_1(mask, :);
plotData_1 = sortrows(plotData_1);
plotData_1 = plotData_1(1:round(length(faceCorresp)*0.99),:);
plotBinData(plotData_1(:,1), plotData_1(:,2), [0, 200], ylim_bin, 1, 'Integral |Face Curvature| in An4, \mum', 'Face Area Difference, \mum^{2}', true, false) 
% print('FItgCurv_An4_FADiff_3', '-dpng','-r300')





