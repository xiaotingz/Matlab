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

% ##### Find Face Correspondence #####
%	--- make a look up list that converts ID_An5 to ID_An4. ---
%	---	if no correspondence, assign NaN ---
lookUp = sortrows(lookUp, 2);  %  somehow the look up table is not really sorted well 
lookUp_5to4 = zeros(max(lookUp(:,2)),2);
%   --- implicitly, grainSizeDiff is for the grains in state_An5 ---
idx = 1;
for i = 1 : max(lookUp(:,2))
    if lookUp(idx,2) == i
        lookUp_5to4(i,1) = lookUp(idx,1);
        lookUp_5to4(i,2) = lookUp(idx,2);
        idx = idx + 1;
    else
        lookUp_5to4(i,1) = NaN;
        lookUp_5to4(i,2) = i;
    end
end
%   --- double check that ID_An5 is implicit correctly  --- 
if sum(lookUp_5to4(:,2) == [1:length(numNeigh_An5)].') == length(numNeigh_An5)
    lookUp_5to4(:,2) = [];
end


% ##### get the faceLabels and their correpondence ##### 
[faces_An4, faces_An5, faceCorresp] = TrackFace(numNeigh_An4, neighList_An4, numNeigh_An5, neighList_An5, lookUp_5to4);


% ##### calc integral face curvature ##### 
%   --- faceCurvs = [integralArea, integralCurvature] ---
%   --- data is the order of faces_An but contain just the valid faces ---
FCurvs_An4 = faceCurvatureForTrack(file_An4, faces_An4);
FCurvs_An5 = faceCurvatureForTrack(file_An5, faces_An5);
diff_tmp = FCurvs_An4(faceCorresp(:,1),:) - FCurvs_An5(faceCorresp(:,2),:);
FA_diff = diff_tmp(:,1);
FCurv_diff = abs(diff_tmp(:,2));
clear diff_tmp


% ##### get the sizeChange of the two grains that defines the face #####
tmp1 = numCells_An4(faces_An4(faceCorresp(:,1),:));
tmp2 = numCells_An5(faces_An5(faceCorresp(:,2),:));
tmp3 = tmp1 - tmp2;
%   --- The face mobility should be associated with delta_G1 - delta_G2, ---
%   --- because a face moves the most when onegrain grow one grain shrink ---
faceMob_dV = abs(tmp3(:,1) - tmp3(:,2));
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
set(0,'defaultAxesLabelFontSize',1.1)
set(0,'defaultAxesFontSize',19)

figure(1)
fig1 = scatter(FCurv_diff,FA_diff,'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
xlabel('Integral Face Curvature Difference')
ylabel('Face Area Difference')
print('FCurvDiff_FADiff', '-dpng','-r300')

figure(2)
scatter(FCurvs_An4(faceCorresp(:,1),2),FA_diff,'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
xlabel('Integral Face Curvature in An4')
ylabel('Face Area Difference')
print('FCurv_An4_FADiff', '-dpng','-r300')

figure(3)
scatter(FCurvs_An5(faceCorresp(:,2),2),FA_diff,'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
xlabel('Integral Face Curvature in An5')
ylabel('Face Area Difference')
print('FCurv_An5_FADiff', '-dpng','-r300')

figure(4)
scatter(FCurvs_An4(faceCorresp(:,1),2),faceMob_dV,'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
xlabel('Integral Face Curvature in An4')
ylabel('The Associated Volume Change')
print('FCurv_An4_Vdiff', '-dpng','-r300')

figure(5)
scatter(FCurv_diff,faceMob_dV,'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
xlabel('Integral Face Curvature Difference')
ylabel('The Associated Volume Change')
print('FCurvDiff_Vdiff', '-dtiff','-r300')

figure(6)
scatter(FCurv_diff,FCurvs_An4(faceCorresp(:,1),1),'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
xlabel('Integral Face Curvature Difference')
ylabel('Face Area in An4')
print('FCurvDiff_FA_An4', '-dpng','-r300')

figure(7)
scatter(FCurv_diff,FCentrs_diff,'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
xlabel('Integral Face Curvature Difference')
ylabel('norm(Face Centorid Difference)')
print('FCurvDiff_CentrDiff', '-dpng','-r300')


% figure(9)
% scatter(misAs_An4(faceCorresp(:,1)),faceMob_dV,'filled');
% xlabel('The Misorientation Angle of Grain Face')
% ylabel('The Associated Volume Change')

% --- binByMisA(property, misAs, gridSize) ---
% --- data_grid = [leftBoundaryOfCurrentBin, avg(property), countInBin] ---
% data_grid = binByMisA(faceMob_dV, misAs_An4(faceCorresp(:,1)), 0.5);
% figure(10)
% scatter(data_grid(:,1), data_grid(:,2),'filled');
% xlabel('Misorientation in An4')
% ylabel('The average dV across face')
