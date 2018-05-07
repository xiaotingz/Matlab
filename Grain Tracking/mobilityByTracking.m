clear 

file_An4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_mesh.dream3d');
file_An5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh.dream3d');
load('lookUpTable_An4_An5.mat')

% centroids_An4 = roundn(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
numCells_An4 = h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/NumElements').';
numNeigh_An4 = double(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
neighList_An4 = double(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
% centroids_An5 = roundn(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
numCells_An5 = h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/NumElements').';
numNeigh_An5 = double(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
neighList_An5 = double(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
% centroids_An4(1,:) = [];
numCells_An4(1) = [];
numNeigh_An4(1) = [];
% centroids_An5(1,:) = [];
numCells_An5(1) = [];
numNeigh_An5(1) = [];


% ##### get the faceLabels and their correpondence ##### 
[faces_An4, faces_An5, faceCorresp] = TrackFace(numNeigh_An4, neighList_An4, numNeigh_An5, neighList_An5, lookUp);

%%
% ##### calc integral face curvature ##### 
%   --- faceCurvs = [integralArea, integralCurvature] ---
%   --- data is the order of faces_An but contain just the valid faces ---
FItgCurvs_An4 = faceCurvatureForTrack(file_An4, faces_An4);
FItgCurvs_An5 = faceCurvatureForTrack(file_An5, faces_An5);
diff_tmp = FItgCurvs_An4(faceCorresp(:,1),:) - FItgCurvs_An5(faceCorresp(:,2),:);
FA_diff = diff_tmp(:,1);
FItgCurv_diff = diff_tmp(:,2);
FAvgCurv_diff = FItgCurvs_An4(faceCorresp(:,1),2)./FItgCurvs_An4(faceCorresp(:,1),1) - FItgCurvs_An5(faceCorresp(:,2),2)./FItgCurvs_An5(faceCorresp(:,2),1);
clear diff_tmp


% ##### get the sizeChange of the two grains that defines the face #####
tmp1 = numCells_An4(faces_An4(faceCorresp(:,1),:));
tmp2 = numCells_An5(faces_An5(faceCorresp(:,2),:));
tmp3 = tmp1 - tmp2;
%   --- The face mobility should be associated with delta_G1 - delta_G2, ---
%   --- because a face moves the most when onegrain grow one grain shrink ---
faceMob_dV = tmp3(:,1) - tmp3(:,2);
clear tmp1 tmp2 tmp3


% ##### get the face coordinates #####
FCentrs_An4 = findFaceCentorids(file_An4, faces_An4);
FCentrs_An5 = findFaceCentorids(file_An5, faces_An5);
FCentrs_diff = FCentrs_An4(faceCorresp(:,1),:) - FCentrs_An5(faceCorresp(:,2),:);
FCentrs_diff = sqrt(sum(FCentrs_diff .* FCentrs_diff, 2));


% ##### get the Rodrigues Vector corresponding to the Faces #####
% [RFvecs_An4, misAs_An4]  = getFaceRFvecs(file_An4, faces_An4);
% [RFvecs_An5, misAs_An5]  = getFaceRFvecs(file_An5, faces_An5);


%%
% --- single plot should use font 1.1 & 19 ---
set(0,'defaultAxesLabelFontSize',1.1)
set(0,'defaultAxesFontSize',14)
% set(0,'defaultAxesFontWeight','bold')

% ----------------------- inputs of plot functions -----------------------
% --- plotScatter(x, y, label_x, label_y)
% --- plotBinData(x, y, xrange, stepsize, xlabel, ylabel, showSTD) 
%         showSTD option is true by default
% ------------------------------------------------------------------------

% ##### Face Area Change #####
figure('rend','painters','pos',[10 10 900 1100])
ylim_bin = [-1e4, 1e4];
subplot(3,2,1)
plotScatter(FItgCurvs_An4(faceCorresp(:,1),2), FA_diff, 'Integral Face Curvature in An4, \mum', 'Face Area Difference, \mum^{2}');
title(['ylim of binning = [', num2str(ylim_bin), '], An4'])
subplot(3,2,3)
plotBinData(FItgCurvs_An4(faceCorresp(:,1),2), FA_diff, [0, 500], ylim_bin, 2, 'Integral Face Curvature in An4, \mum', 'Face Area Difference, \mum^{2}', true, false) 
subplot(3,2,5)
plotBinData(FItgCurvs_An4(faceCorresp(:,1),2), FA_diff, [0, 200], ylim_bin, 2, 'Integral Face Curvature in An4, \mum', 'Face Area Difference, \mum^{2}', true, true)
subplot(3,2,2)
plotScatter(FItgCurvs_An5(faceCorresp(:,2),2), FA_diff, 'Integral Face Curvature in An5, \mum', 'Face Area Difference, \mum^{2}');
title(['ylim of binning = [', num2str(ylim_bin), '], An5'])
subplot(3,2,4)
plotBinData(FItgCurvs_An5(faceCorresp(:,2),2), FA_diff, [0, 500], ylim_bin, 2, 'Integral Face Curvature in An5, \mum', 'Face Area Difference, \mum^{2}', true, false)
subplot(3,2,6)
plotBinData(FItgCurvs_An5(faceCorresp(:,2),2), FA_diff, [0, 200], ylim_bin, 2, 'Integral Face Curvature in An5, \mum', 'Face Area Difference, \mum^{2}', true, true)
print('FItgCurv_An5_FADiff', '-dpng','-r300')
%%
figure('rend','painters','pos',[10 10 400 1100])
ylim_bin = [-1e4, 1e4];
subplot(3,1,1)
plotScatter(FItgCurv_diff, FA_diff, 'Integral Face Curvature Difference, \mum', 'Face Area Difference, \mum^{2}');
title(['ylim of binning = [', num2str(ylim_bin), ']'])
subplot(3,1,2)
plotBinData(FItgCurv_diff, FA_diff, [-500, 500], ylim_bin, 2, 'Integral Face Curvature Difference, \mum', 'Face Area Difference, \mum^{2}', true, true )
subplot(3,1,3)
plotBinData(FItgCurv_diff, FA_diff, [-100, 100], ylim_bin, 2, 'Integral Face Curvature Difference, \mum', 'Face Area Difference, \mum^{2}', true, true)
print('FItgCurvDiff_FADiff', '-dpng','-r300')

%%
figure('rend','painters','pos',[10 10 1300 1100])
ylim_bin = [-1e4, 1e4];
subplot(2,2,1)
plotScatter(FAvgCurv_diff, FA_diff, 'Average Face Curvature Difference, \mum^{-1}', 'Face Area Difference, \mum^{2}');
title(['ylim of binning = [', num2str(ylim_bin), ']'])
xlim([-50, 50])
subplot(2,2,2)
plotBinData(FAvgCurv_diff, FA_diff, [-20, 20], ylim_bin, 0.5, 'Average Face Curvature Difference, \mum^{-1}', 'Face Area Difference, \mum^{2}', true, true)
subplot(2,2,3)
plotBinData(FAvgCurv_diff, FA_diff, [-6, 6], ylim_bin, 0.05, 'Average Face Curvature Difference, \mum^{-1}', 'Face Area Difference, \mum^{2}', true, false)
subplot(2,2,4)
plotBinData(FAvgCurv_diff, FA_diff, [-0.5, 0.5], ylim_bin, 0.005, 'Average Face Curvature Difference, \mum^{-1}', 'Face Area Difference, \mum^{2}', true, false)
% ylim([-1000, 1000])
print('FAvgCurvDiff_FADiff', '-dpng','-r300')
%%
% ##### Associated Volume Change #####
figure('rend','painters','pos',[10 10 1300 1100])
ylim_bin = [0, 14e4];
subplot(2,2,1)
plotScatter(FItgCurvs_An4(faceCorresp(:,1),2), abs(faceMob_dV), 'Integral Face Curvature in An4, \mum', 'The Associated Volume Change, \mum^{3}');
title(['ylim of binning = [', num2str(ylim_bin), ']'])
subplot(2,2,2)
plotBinData(FItgCurvs_An4(faceCorresp(:,1),2), abs(faceMob_dV), [0, 1000], ylim_bin, 10, 'Integral Face Curvature in An4, \mum', 'The Associated Volume Change, \mum^{3}', true, false)
subplot(2,2,3)
plotBinData(FItgCurvs_An4(faceCorresp(:,1),2), abs(faceMob_dV), [0, 150], ylim_bin, 5, 'Integral Face Curvature in An4, \mum', 'The Associated Volume Change, \mum^{3}', true, false)
print('FItgCurv_An4_Vdiff', '-dpng','-r300')
subplot(2,2,4)
plotBinData(FItgCurvs_An4(faceCorresp(:,1),2), abs(faceMob_dV), [0, 150], ylim_bin, 5, 'Integral Face Curvature in An4, \mum', 'The Associated Volume Change, \mum^{3}', true, true)
print('FItgCurv_An4_Vdiff', '-dpng','-r300')

%%
% ##### Face Integral Curvature Change #####
figure('rend','painters','pos',[10 10 1300 1100])
ylim_bin = [-1e4, 1e4];
subplot(2,2,1)
plotScatter(FItgCurvs_An4(faceCorresp(:,1),2), FItgCurv_diff, 'Integral Face Curvature in An4, \mum', 'Integral Face Curvature Difference, \mum^{2}');
title(['ylim of binning = [', num2str(ylim_bin), ']'])
subplot(2,2,2)
plotBinData(FItgCurvs_An4(faceCorresp(:,1),2), FItgCurv_diff, [0, 300], ylim_bin, 5, 'Integral Face Curvature in An4, \mum', 'Integral Face Curvature Difference, \mum^{2}', false)
subplot(2,2,3)
plotScatter(FItgCurvs_An5(faceCorresp(:,2),2), FItgCurv_diff, 'Integral Face Curvature in An5, \mum', 'Integral Face Curvature Difference, \mum^{2}');
title(['ylim of binning = [', num2str(ylim_bin), ']'])
subplot(2,2,4)
plotBinData(FItgCurvs_An5(faceCorresp(:,2),2), FItgCurv_diff, [0, 300], ylim_bin, 5, 'Integral Face Curvature in An5, \mum', 'Integral Face Curvature Difference, \mum^{2}', false)
print('FItgCurv_FItgCurvDiff', '-dpng','-r300')

% ##### Face Centroid Position Change #####
%%
figure('rend','painters','pos',[10 10 1300 1100])
ylim_bin = [0, 1e4];
subplot(2,2,1)
plotScatter(abs(FItgCurv_diff), FCentrs_diff, 'abs(Integral Face Curvature Difference), \mum', 'norm(Face Centorid Difference), \mum');
title(['ylim of binning = [', num2str(ylim_bin), ']'])
subplot(2,2,2)
plotBinData(abs(FItgCurv_diff), FCentrs_diff, [0, 500], ylim_bin, 5, 'abs(Integral Face Curvature Difference), \mum', 'norm(Face Centorid Difference), \mum')
subplot(2,2,3)
plotBinData(abs(FItgCurv_diff), FCentrs_diff, [0, 200], ylim_bin, 2, 'abs(Integral Face Curvature Difference), \mum', 'norm(Face Centorid Difference), \mum')
subplot(2,2,4)
plotBinData(abs(FItgCurv_diff), FCentrs_diff, [0, 200], ylim_bin, 2, 'abs(Integral Face Curvature Difference), \mum', 'norm(Face Centorid Difference), \mum', false)
print('FItgCurvDiff_CentrDiff', '-dpng','-r300')
%%

% figure(9)
% scatter(misAs_An4(faceCorresp(:,1)),faceMob_dV,'filled');
% xlabel('The Misorientation Angle of Grain Face')
% ylabel('The Associated Volume Change')

% --- binByMisA(property, misAs, gridSize) ---
% --- data_grid = [leftBoundaryOfCurrentBin, avg(property), countInBin] ---
% ### SHOULD TRY DO IT BY binDATA ###
% data_grid = binByMisA(faceMob_dV, misAs_An4(faceCorresp(:,1)), 0.5);
% figure(10)
% scatter(data_grid(:,1), data_grid(:,2),'filled');
% xlabel('Misorientation in An4')
% ylabel('The average dV across face')
