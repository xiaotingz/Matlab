clear 

file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
load('look_up_table_an4_an5.mat')

% centroids_An4 = roundn(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
num_cells_an4 = h5read(file_an4,'/DataContainers/ImageDataContainer/CellFeatureData/NumElements').';
% centroids_An5 = roundn(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
num_cells_an5 = h5read(file_an5,'/DataContainers/ImageDataContainer/CellFeatureData/NumElements').';
% centroids_An4(1,:) = [];
num_cells_an4(1) = [];
% centroids_An5(1,:) = [];
num_cells_an5(1) = [];

% ##### get the faceLabels and their correpondence ##### 
[faces_an4, faces_an5, face_corresp] = TrackFace(file_an4, file_an5, look_up_table, true);

%%
% ##### calc integral |face curvature| ##### 
%   --- faceCurvs = [integralArea, integralCurvature] ---
%   --- data is the order of faces_An but contain just the valid faces ---
face_integ_curv_an4 = FaceCurvatureForTrack(file_an4, faces_an4);
face_integ_curv_an5 = FaceCurvatureForTrack(file_an5, faces_an5);
diff_tmp = face_integ_curv_an5(face_corresp(:,2),:) - face_integ_curv_an4(face_corresp(:,1),:);
face_area_diff = diff_tmp(:,1);
face_integ_curv_diff = diff_tmp(:,2);
diff = face_integ_curv_an5(face_corresp(:,2),2)./face_integ_curv_an5(face_corresp(:,2),1) - face_integ_curv_an4(face_corresp(:,1),2)./face_integ_curv_an4(face_corresp(:,1),1);
clear diff_tmp


% ##### get the sizeChange of the two grains that defines the face #####
tracked_facelabel_an4 = faces_an4(face_corresp(:,1),:);
tracked_facelabel_an5 = faces_an5(face_corresp(:,2),:);
map_5from4 = containers.Map(look_up_table(:,1), look_up_table(:,2));
for i = 1:length(tracked_facelabel_an4)
    if tracked_facelabel_an5(i,1) ~= map_5from4(tracked_facelabel_an4(i,1))
        tracked_facelabel_an5(i,:) = flip(tracked_facelabel_an5(i,:));
    end
end
tmp1 = num_cells_an4(tracked_facelabel_an4);
tmp2 = num_cells_an5(tracked_facelabel_an5);
tmp3 = tmp2 - tmp1;
%   --- The face mobility should be associated with delta_G1 - delta_G2, ---
%   --- because a face moves the most when onegrain grow one grain shrink ---
% cellVolume = 2.8125*2.8125*4;
faceMob_dV = abs(tmp3(:,1) - tmp3(:,2));
clear tmp1 tmp2 tmp3


% ##### get the face coordinates #####
face_centroid_an4 = findFaceCentorids(file_an4, faces_an4);
face_centroid_an5 = findFaceCentorids(file_an5, faces_an5);
face_centroid_diff = face_centroid_an5(face_corresp(:,2),:) - face_centroid_an4(face_corresp(:,1),:);
face_centroid_diff = sqrt(sum(face_centroid_diff .* face_centroid_diff, 2));


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

% ##### Face Area v.s. Face Integral Curvature #####
figure('rend','painters','pos',[10 10 900 1100])
ylim_bin = [-1e4, 1e4];
subplot(3,2,1)
plotScatter(face_integ_curv_an4(face_corresp(:,1),1), face_integ_curv_an4(face_corresp(:,1),2), 'Face Area in An4, \mum^{2}', 'Integral |Face Curvature| in An4, \mum');
title(['ylim of binning = [', num2str(ylim_bin), '], An4'])
subplot(3,2,3)
plotBinData(face_integ_curv_an4(face_corresp(:,1),1), face_integ_curv_an4(face_corresp(:,1),2), [0, 10000], ylim_bin, 100, 'Face Area in An4, \mum^{2}', 'Integral |Face Curvature| in An4, \mum', true, false) 
subplot(3,2,5)
plotBinData(face_integ_curv_an4(face_corresp(:,1),1), face_integ_curv_an4(face_corresp(:,1),2), [0, 1000], ylim_bin, 50, 'Face Area in An4, \mum^{2}', 'Integral |Face Curvature| in An4, \mum', true, true)
subplot(3,2,2)
plotScatter(face_integ_curv_an5(face_corresp(:,2),1), face_integ_curv_an5(face_corresp(:,2),2), 'Face Area in An5, \mum^{2}', 'Integral |Face Curvature| in An5, \mum');
title(['ylim of binning = [', num2str(ylim_bin), '], An5'])
subplot(3,2,4)
plotBinData(face_integ_curv_an5(face_corresp(:,2),1), face_integ_curv_an5(face_corresp(:,2),2), [0, 10000], ylim_bin, 100, 'Face Area in An5, \mum^{2}', 'Integral |Face Curvature| in An5, \mum', true, false)
subplot(3,2,6)
plotBinData(face_integ_curv_an5(face_corresp(:,2),1), face_integ_curv_an5(face_corresp(:,2),2), [0, 1000], ylim_bin, 50, 'Face Area in An5, \mum^{2}', 'Integral |Face Curvature| in An5, \mum', true, true)
print('FArea_FIntCurv', '-dpng','-r300')

%% ##### Face Integral Curvature v.s. Face Area Change #####
figure('rend','painters','pos',[10 10 900 1100])
ylim_bin = [-1e4, 1e4];
subplot(3,2,1) 
plotScatter(face_integ_curv_an4(face_corresp(:,1),2), face_area_diff, 'Integral |Face Curvature| in An4, \mum', 'Face Area Difference, \mum^{2}');
title(['ylim of binning = [', num2str(ylim_bin), '], An4'])
subplot(3,2,3)
plotBinData(face_integ_curv_an4(face_corresp(:,1),2), face_area_diff, [0, 500], ylim_bin, 2, 'Integral |Face Curvature| in An4, \mum', 'Face Area Difference, \mum^{2}', true, false) 
subplot(3,2,5)
plotBinData(face_integ_curv_an4(face_corresp(:,1),2), face_area_diff, [0, 200], ylim_bin, 2, 'Integral |Face Curvature| in An4, \mum', 'Face Area Difference, \mum^{2}', true, true)
subplot(3,2,2)
plotScatter(face_integ_curv_an5(face_corresp(:,2),2), face_area_diff, 'Integral |Face Curvature| in An5, \mum', 'Face Area Difference, \mum^{2}');
title(['ylim of binning = [', num2str(ylim_bin), '], An5'])
subplot(3,2,4)
plotBinData(face_integ_curv_an5(face_corresp(:,2),2), face_area_diff, [0, 500], ylim_bin, 2, 'Integral |Face Curvature| in An5, \mum', 'Face Area Difference, \mum^{2}', true, false)
subplot(3,2,6)
plotBinData(face_integ_curv_an5(face_corresp(:,2),2), face_area_diff, [0, 200], ylim_bin, 2, 'Integral |Face Curvature| in An5, \mum', 'Face Area Difference, \mum^{2}', true, true)
print('FItgCurv_FADiff', '-dpng','-r300')

%%
figure('rend','painters','pos',[10 10 400 1100])
ylim_bin = [-1e4, 1e4];
subplot(3,1,1)
plotScatter(face_integ_curv_diff, face_area_diff, 'Integral |Face Curvature| Difference, \mum', 'Face Area Difference, \mum^{2}');
title(['ylim of binning = [', num2str(ylim_bin), ']'])
subplot(3,1,2)
plotBinData(face_integ_curv_diff, face_area_diff, [-500, 500], ylim_bin, 5, 'Integral |Face Curvature| Difference, \mum', 'Face Area Difference, \mum^{2}', true, true )
subplot(3,1,3)
plotBinData(face_integ_curv_diff, face_area_diff, [-100, 100], ylim_bin, 5, 'Integral |Face Curvature| Difference, \mum', 'Face Area Difference, \mum^{2}', true, true)
print('FItgCurvDiff_FADiff', '-dpng','-r300')

%% ###### Average Face Curvature Difference v.s. Face Area Difference #####
figure('rend','painters','pos',[10 10 1300 1100])
ylim_bin = [-1e4, 1e4];
subplot(2,2,1)
plotScatter(FAvgCurv_diff, face_area_diff, 'Average Face Curvature Difference, \mum^{-1}', 'Face Area Difference, \mum^{2}');
title(['ylim of binning = [', num2str(ylim_bin), ']'])
% xlim([-50, 50])
subplot(2,2,2)
plotBinData(FAvgCurv_diff, face_area_diff, [-20, 20], ylim_bin, 0.5, 'Average Face Curvature Difference, \mum^{-1}', 'Face Area Difference, \mum^{2}', true, true)
subplot(2,2,3)
plotBinData(FAvgCurv_diff, face_area_diff, [-6, 6], ylim_bin, 0.05, 'Average Face Curvature Difference, \mum^{-1}', 'Face Area Difference, \mum^{2}', true, false)
subplot(2,2,4)
plotBinData(FAvgCurv_diff, face_area_diff, [-0.5, 0.5], ylim_bin, 0.01, 'Average Face Curvature Difference, \mum^{-1}', 'Face Area Difference, \mum^{2}', true, false)
% ylim([-1000, 1000])
print('FAvgCurvDiff_FADiff', '-dpng','-r300')

%% ##### Integral Face Curvature in An4 v.s. Integral Face Curvature Difference #####
figure('rend','painters','pos',[10 10 1300 1100])
ylim_bin = [-1e4, 1e4];
subplot(2,2,1)
plotScatter(face_integ_curv_an4(face_corresp(:,1),2), face_integ_curv_diff, 'Integral |Face Curvature| in An4, \mum', 'Integral |Face Curvature| Difference, \mum');
title(['ylim of binning = [', num2str(ylim_bin), ']'])
subplot(2,2,2)
plotBinData(face_integ_curv_an4(face_corresp(:,1),2), face_integ_curv_diff, [0, 300], ylim_bin, 5, 'Integral |Face Curvature| in An4, \mum', 'Integral |Face Curvature| Difference, \mum', false)
subplot(2,2,3)
plotScatter(face_integ_curv_an5(face_corresp(:,2),2), face_integ_curv_diff, 'Integral |Face Curvature| in An5, \mum', 'Integral |Face Curvature| Difference, \mum');
title(['ylim of binning = [', num2str(ylim_bin), ']'])
subplot(2,2,4)
plotBinData(face_integ_curv_an5(face_corresp(:,2),2), face_integ_curv_diff, [0, 300], ylim_bin, 5, 'Integral |Face Curvature| in An5, \mum', 'Integral |Face Curvature| Difference, \mum', false)
print('FItgCurv_FItgCurvDiff', '-dpng','-r300')

%% ##### Integral Face Curvature Difference v.s. Associated Volume Difference  #####
figure('rend','painters','pos',[10 10 1300 1100])
ylim_bin = [0, 14e4];
subplot(2,2,1)
plotScatter(face_integ_curv_an4(face_corresp(:,1),2), faceMob_dV, 'Integral |Face Curvature| in An4, \mum', 'The Associated Volume Change, \mum^{3}');
title(['ylim of binning = [', num2str(ylim_bin), ']'])
subplot(2,2,2)
plotBinData(face_integ_curv_an4(face_corresp(:,1),2), faceMob_dV, [0, 1000], ylim_bin, 5, 'Integral |Face Curvature| in An4, \mum', 'The Associated Volume Change, \mum^{3}', true, false)
subplot(2,2,3)
plotBinData(face_integ_curv_an4(face_corresp(:,1),2), faceMob_dV, [0, 200], ylim_bin, 10, 'Integral |Face Curvature| in An4, \mum', 'The Associated Volume Change, \mum^{3}', true, false)
print('FItgCurv_An4_Vdiff', '-dpng','-r300')
subplot(2,2,4)
plotBinData(face_integ_curv_an4(face_corresp(:,1),2), faceMob_dV, [0, 150], ylim_bin, 5, 'Integral |Face Curvature| in An4, \mum', 'The Associated Volume Change, \mum^{3}', true, true)
print('FItgCurv_An4_Vdiff', '-dpng','-r300')

%% ##### |Integral Face Curvature Difference| v.s. Face Centroid Position Difference #####
figure('rend','painters','pos',[10 10 1300 1100])
ylim_bin = [0, 1e4];
subplot(2,2,1)
plotScatter(abs(face_integ_curv_diff), face_centroid_diff, 'abs(Integral |Face Curvature| Difference), \mum', 'norm(Face Centorid Difference), \mum');
title(['ylim of binning = [', num2str(ylim_bin), ']'])
subplot(2,2,2)
plotBinData(abs(face_integ_curv_diff), face_centroid_diff, [0, 500], ylim_bin, 5, 'abs(Integral |Face Curvature| Difference), \mum', 'norm(Face Centorid Difference), \mum',true,true)
subplot(2,2,3)
plotBinData(abs(face_integ_curv_diff), face_centroid_diff, [0, 200], ylim_bin, 2, 'abs(Integral |Face Curvature| Difference), \mum', 'norm(Face Centorid Difference), \mum',true, true)
subplot(2,2,4)
plotBinData(abs(face_integ_curv_diff), face_centroid_diff, [0, 200], ylim_bin, 2, 'abs(Integral |Face Curvature| Difference), \mum', 'norm(Face Centorid Difference), \mum', true, false)
% ylim([450,750])
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
