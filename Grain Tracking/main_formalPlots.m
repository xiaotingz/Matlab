set(0,'defaultAxesLabelFontSize',1.2)
set(0,'defaultAxesFontSize',17)

ylim_bin = [-100000, 100000];

%% ##### Integral Face Curvature in An4 v.s. Face Area Difference #####
plotData_1 = [FItgCurvs_An4(faceCorresp(:,1),2), FA_diff];
plotData_1 = filterExtremeData(plotData_1, 0.99);
plotBinData(plotData_1(:,1), plotData_1(:,2), [0, 420], ylim_bin, 5, 'Integral |Face Curvature| in An4, \mum', 'Face Area Difference, \mum^{2}', true, false) 
print('FItgCurvAn4_FADiff_3', '-dpng','-r300')

%% ##### Average Face Curvature Difference v.s. Face Area Difference #####
plotData_2 = [FAvgCurv_diff, FA_diff];
plotData_2 = filterExtremeData(plotData_2, 1);
plotBinData(plotData_2(:,1), plotData_2(:,2),[-0.6, 0.6], ylim_bin, 0.02, 'Average Face Curvature Difference, \mum^{-1}', 'Face Area Difference, \mum^{2}', true, false)
print('FAvgCurvDiff_FADiff_1', '-dpng','-r300')

%% ##### Integral Face Curvature in An4 v.s. Integral Face Curvature Difference #####
plotData_3 = [FItgCurvs_An4(faceCorresp(:,1),2), FItgCurv_diff];
plotData_3 = filterExtremeData(plotData_3, 1);
plotBinData(plotData_3(:,1), plotData_3(:,2), [0, 300], ylim_bin, 5, 'Integral |Face Curvature| in An4, \mum', 'Integral |Face Curvature| Difference, \mum', true, true)
print('FItgCurvsAn4__FItgCurvDiff', '-dpng','-r300')

%% ##### Integral Face Curvature Difference v.s. Associated Volume Difference  #####
plotData_4 = [FItgCurvs_An4(faceCorresp(:,1),2), abs(faceMob_dV)];
plotData_4 = filterExtremeData(plotData_4, 0.99);
plotBinData(plotData_4(:,1), plotData_4(:,2), [0, 200], ylim_bin, 10, 'Integral |Face Curvature| in An4, \mum', 'The Associated Volume Change, \mum^{3}', false, true)
% ylim([1200,2200])
print('FItgCurvsAn4_|faceMob_dV|_2', '-dpng','-r300')

%% ##### |Integral Face Curvature Difference| v.s. Face Centroid Position Difference  #####
plotData_5 = [abs(FItgCurv_diff), FCentrs_diff];
plotData_5 = filterExtremeData(plotData_5, 0.99);
plotBinData(plotData_5(:,1), plotData_5(:,2), [0, 200], ylim_bin, 5, 'abs(Integral |Face Curvature| Difference), \mum', 'norm(Face Centorid Difference), \mum', false, false)
print('|FItgCurv_diff|_FCentrs_diff_2', '-dpng','-r300')



%%
C_diff = numEdges_An5(faceCorresp(:,2)) - numEdges_An4(faceCorresp(:,1));
C_An4 = numEdges_An4(faceCorresp(:,1));
C_An5 = numEdges_An5(faceCorresp(:,2));
%% ##### #(Face Corners) Difference v.s. Integral Face Curvature Difference #####
plotData_6 = [C_diff, FA_diff];
plotData_6 = filterExtremeData(plotData_6, 1);
plotBinData(plotData_6(:,1), plotData_6(:,2), [-30, 20], ylim_bin, 1, '#(Face Corners) Difference', 'Face Area Difference, \mum^{2}', true, true) 
print('CDiff_FADiff_1', '-dpng','-r300')

%% ##### #(Face Corners) Difference v.s. Associated Volume Difference #####
plotData_7 = [C_diff, abs(faceMob_dV)];
plotData_7 = filterExtremeData(plotData_7, 0.99);
plotBinData(plotData_7(:,1), plotData_7(:,2), [-10, 10], ylim_bin, 1, '#(Face Corners) Difference', 'The Associated Volume Change, \mum^{3}', false, true)
print('CDiff_|faceMob_dV|_2', '-dpng','-r300')

%% ##### #(Face Corners) Difference v.s. Face Centroid Position Difference #####
plotData_8 = [C_diff, FCentrs_diff];
plotData_8 = filterExtremeData(plotData_8, 0.99);
plotBinData(plotData_8(:,1), plotData_8(:,2), [-10, 10], ylim_bin, 1, '#(Face Corners) Difference', 'norm(Face Centorid Difference), \mum', false, true)
print('CDiff_FCentrs_diff_2', '-dpng','-r300')

