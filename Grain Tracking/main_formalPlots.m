set(0,'defaultAxesLabelFontSize',1.2)
set(0,'defaultAxesFontSize',17)

%% ##### The Complete Tracked Grains #####
SG_An4 = h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
SG_An5 = h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
numNeigh_An4 = double(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
numNeigh_An5 = double(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
D_An4 = double(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters')).';
D_An5 = double(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters')).';
SG_An4(1) = [];
SG_An5(1) = [];
numNeigh_An4(1) = [];
numNeigh_An5(1) = [];
D_An4(1) = [];
D_An5(1) = [];

mask_lookUp = logical([SG_An4(lookUp(:,1)), SG_An5(lookUp(:,2))]);
lookUp_CG = lookUp(all(~mask_lookUp, 2), :);

% ##### Log Normal Size Distribution #####
D_An4_TG = D_An4(lookUp(:,1));
D_An5_TG = D_An5(lookUp(:,2));
D_An4_TCG = D_An4(lookUp_CG(:,1));
D_An5_TCG = D_An5(lookUp_CG(:,2));
stepsize = 0.2;
grid_D_A4TG = BinData(log(D_An4_TG/(sum(D_An4_TG)/length(D_An4_TG))), stepsize);
grid_D_A5TG = BinData(log(D_An5_TG/(sum(D_An5_TG)/length(D_An5_TG))), stepsize);
grid_D_A4TCG = BinData(log(D_An4_TCG/(sum(D_An4_TCG)/length(D_An4_TCG))), stepsize);
grid_D_A5TCG = BinData(log(D_An5_TCG/(sum(D_An5_TCG)/length(D_An5_TCG))), stepsize);
figure(1)
plot(grid_D_A4TG(:,1), grid_D_A4TG(:,2),'-o','Color','k','MarkerFaceColor', 'k','LineWidth',1.5);
hold on
plot(grid_D_A5TG(:,1), grid_D_A5TG(:,2),'-^','Color','r','MarkerFaceColor', 'r','LineWidth',1.5);
plot(grid_D_A4TCG(:,1), grid_D_A4TCG(:,2),'--o','Color','k','LineWidth',1.5);
plot(grid_D_A5TCG(:,1), grid_D_A5TCG(:,2),'--^','Color','r','LineWidth',1.5);
xlim([-1.5, 2.5]);
legend('Anneal 4, all tracked grains', 'Anneal 5, all tracked grains', 'Anneal 4, tracked inner grains', 'Anneal 5, tracked inner grains')
xlabel('log(D/<D>)');
ylabel('Frequency');
print('Size log normal distri', '-dpng','-r300')

% ##### #Faces Distribution #####
numNeigh_An4_TG = numNeigh_An4(lookUp(:,1));
numNeigh_An5_TG = numNeigh_An5(lookUp(:,2));
numNeigh_An4_TCG = numNeigh_An4(lookUp_CG(:,1));
numNeigh_An5_TCG = numNeigh_An5(lookUp_CG(:,2));
stepsize = 1;
grid_numNeigh_A4TG = BinData(numNeigh_An4_TG, stepsize);
grid_numNeigh_A5TG = BinData(numNeigh_An5_TG, stepsize);
grid_numNeigh_A4TCG = BinData(numNeigh_An4_TCG, stepsize);
grid_numNeigh_A5TCG = BinData(numNeigh_An5_TCG, stepsize);
figure()
% plot(grid_numNeigh_A4TG(:,1), grid_numNeigh_A4TG(:,2),'-o','Color','k','MarkerFaceColor', 'k','LineWidth',1.5);
% hold on
% plot(grid_numNeigh_A5TG(:,1), grid_numNeigh_A5TG(:,2),'-^','Color','r','MarkerFaceColor', 'r','LineWidth',1.5);
% plot(grid_numNeigh_A4TCG(:,1), grid_numNeigh_A4TCG(:,2),'--o','Color','k','LineWidth',1.5);
% plot(grid_numNeigh_A5TCG(:,1), grid_numNeigh_A5TCG(:,2),'--^','Color','r','LineWidth',1.5);
plot(grid_numNeigh_A4TG(:,1), grid_numNeigh_A4TG(:,2),'Color','k','LineWidth',1.5);
hold on
plot(grid_numNeigh_A5TG(:,1), grid_numNeigh_A5TG(:,2),'Color','r','LineWidth',1.5);
plot(grid_numNeigh_A4TCG(:,1), grid_numNeigh_A4TCG(:,2),'--','Color','k','LineWidth',1.5);
plot(grid_numNeigh_A5TCG(:,1), grid_numNeigh_A5TCG(:,2),'--','Color','r','LineWidth',1.5);

legend('Anneal 4, all tracked grains', 'Anneal 5, all tracked grains', 'Anneal 4, tracked inner grains', 'Anneal 5, tracked inner grains')
xlabel('F');
ylabel('Frequency');
print('numFaces distri', '-dpng','-r300')



%%
% ####################################################################################################


ylim_bin = [-100000, 100000];

% -----------------------------------------------------------
% filteredData = filterExtremeData(inputData, keepRatio)
%       the final keepRatio = 1 - (1 - keepRatio)*2
% -----------------------------------------------------------
%% ##### Integral Face Curvature in An4 v.s. Face Area Difference #####
plotData_1 = [FItgCurvs_An4(faceCorresp(:,1),2), FA_diff];
plotData_1 = filterExtremeData(plotData_1, 1);
plotBinData(plotData_1(:,1), plotData_1(:,2), [0, 300], ylim_bin, 10, '$|\mathcal{H}|^{F}_{An4}, \mu m$', '$\Delta \mathcal{A}^{F}, \mu m^{2}$', true, false) 
print('FItgCurvAn4_FADiff_allFaces_2', '-dpng','-r300')

%% ##### Average Face Curvature Difference v.s. Face Area Difference #####
plotData_2 = [FAvgCurv_diff, FA_diff];
plotData_2 = filterExtremeData(plotData_2, 1);
plotBinData(plotData_2(:,1), plotData_2(:,2),[-1, 1], ylim_bin, 0.04, '$\Delta <|H|^{F}>, \mu m^{-1}$', '$\Delta \mathcal{A}^{F}, \mu m^{2}$', true, true)
print('FAvgCurvDiff_FADiff_1', '-dpng','-r300')

%% ##### Integral Face Curvature in An4 v.s. Integral Face Curvature Difference #####
plotData_3 = [FItgCurvs_An4(faceCorresp(:,1),2), FItgCurv_diff];
plotData_3 = filterExtremeData(plotData_3, 1);
plotBinData(plotData_3(:,1), plotData_3(:,2), [0, 500], ylim_bin, 10, '$|\mathcal{H}|^{F}_{An4}, \mu m$', '$\Delta |\mathcal{H}|^{F}, \mu m$', true, true)
print('FItgCurvsAn4__FItgCurvDiff_1', '-dpng','-r300')

%% ##### Integral Face Curvature v.s. Associated Volume Difference  #####
plotData_4 = [FItgCurvs_An4(faceCorresp(:,1),2), abs(faceMob_dV)];
plotData_4 = filterExtremeData(plotData_4, 1);
plotBinData(plotData_4(:,1), plotData_4(:,2), [0, 300], [-8e4, 8e4], 10, '$|\mathcal{H}|^{F}_{An4}$, $\mu m$', '$\Delta V^{F}$, \#voxels', true, false)
% ylim([1200,2200])
print('FItgCurvsAn4_|faceMob_dV|_2', '-dpng','-r300')

%% ##### |Integral Face Curvature Difference| v.s. Face Centroid Position Difference  #####
plotData_5 = [FItgCurv_diff, FCentrs_diff];
plotData_5 = filterExtremeData(plotData_5, 1);
plotBinData(plotData_5(:,1), plotData_5(:,2), [-300, 300], ylim_bin, 10, '$\Delta |\mathcal{H}|^{F}, \mu m$', '$||\Delta Centroid^{F}||_{2}, \mu m$', true, true)
print('|FItgCurv_diff|_FCentrs_diff_1', '-dpng','-r300')



%%
C_diff = C_An5(faceCorresp(:,2)) - C_An4(faceCorresp(:,1));
C_An4 = C_An4(faceCorresp(:,1));
C_An5 = C_An5(faceCorresp(:,2));
%% ##### #(Face Corners) Difference v.s. Face Area Difference #####
plotData_6 = [C_diff, FA_diff];
plotData_6 = filterExtremeData(plotData_6, 1);
plotBinData(plotData_6(:,1), plotData_6(:,2), [-10, 10], ylim_bin, 1, '$\Delta C^{F}$', '$\Delta \mathcal{A}^{F}, \mu m^{2}$', true, false) 
print('CDiff_FADiff_2', '-dpng','-r300')

%% ##### #(Face Corners) Difference v.s. Associated Volume Difference #####
plotData_7 = [C_diff, abs(faceMob_dV)];
plotData_7 = filterExtremeData(plotData_7, 1);
plotBinData(plotData_7(:,1), plotData_7(:,2), [-25, 20], ylim_bin, 1, '$\Delta C^{F}$', '$\Delta V^{F}$, \#voxels', true, true)
print('CDiff_|faceMob_dV|_1', '-dpng','-r300')


%% #####  #(Face Corners) in An4 v.s. Face Centroid Position Difference #####
plotData_8 = [C_An4, FA_diff];
plotData_8 = filterExtremeData(plotData_8, 1);
plotBinData(plotData_8(:,1), plotData_8(:,2), [-1, 35], ylim_bin, 1, '$C^{F}_{An4}$', '$\Delta \mathcal{A}^{F}, \mu m^{2}$', true, true)
print('CAn4_FAdiff_1', '-dpng','-r300')



%% ##### #(Face Corners) Difference v.s. Face Centroid Position Difference #####
plotData_8 = [C_diff, FCentrs_diff];
plotData_8 = filterExtremeData(plotData_8, 1);
plotBinData(plotData_8(:,1), plotData_8(:,2), [-10, 10], ylim_bin, 1, '$\Delta C^{F}$', '$||\Delta Centroid^{F}||_{2}, \mu m$', true, false)
print('CDiff_FCentrs_diff_2', '-dpng','-r300')

