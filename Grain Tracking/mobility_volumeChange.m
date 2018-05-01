clear 

file_An4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_mesh.dream3d');
file_An5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh.dream3d');
load('lookUpTable_An4_An5.mat')

% --- feature data ---
% centroids_An4 = roundn(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
numCells_An4 = h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/NumElements').';
numNeigh_An4 = double(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
neighList_An4 = double(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
% centroids_An5 = roundn(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
numCells_An5 = h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/NumElements').';
numNeigh_An5 = double(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
neighList_An5 = double(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
% --- mesh geometries ---
% triNodes_An4 = double(h5read(file_An4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'));
% NCoords_An4 = double(h5read(file_An4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
% triNodes_An5 = double(h5read(file_An5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'));
% NCoords_An5 = double(h5read(file_An5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
% centroids_An4(1,:) = [];
numCells_An4(1) = [];
numNeigh_An4(1) = [];
% centroids_An5(1,:) = [];
numCells_An5(1) = [];
numNeigh_An5(1) = [];

% --- Find Face Correspondence --- 
%	### make a look up list that converts ID_An5 to ID_An4. ###
%	###	if no correspondence, assign NaN ### 
lookUp = sortrows(lookUp, 2);  %  somehow the look up table is not really sorted well 
lookUp_5to4 = zeros(max(lookUp(:,2)),2);
%   ### implicitly, grainSizeDiff is for the grains in state_An5 ###
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
%   ### double check that ID_An5 is implicit correctly  ### 
if sum(lookUp_5to4(:,2) == [1:length(numNeigh_An5)].') == length(numNeigh_An5)
    lookUp_5to4(:,2) = [];
end

% -- get the faceLabels and their correpondence -- 
[faces_An4, faces_An5, faceCorresp] = TrackFace(numNeigh_An4, neighList_An4, numNeigh_An5, neighList_An5, lookUp_5to4);

%%
% -- calc integral face curvature -- 
%   ### faceCurvs = [integralArea, integralCurvature] ###
%   ### data is the order of faces_An but contain just the valid faces ###
faceCurves_An4 = faceCurvatureForTrack(file_An4, faces_An4);
faceCurves_An5 = faceCurvatureForTrack(file_An5, faces_An5);
diff_tmp = faceCurves_An4(faceCorresp(:,1),:) - faceCurves_An5(faceCorresp(:,2),:);
face_areaDiff = diff_tmp(:,1);
face_curvDiff = abs(diff_tmp(:,2));
clear diff_tmp

% --- get the sizeChange of the two grains that defines the face ---
tmp1 = numCells_An4(faces_An4(faceCorresp(:,1),:));
tmp2 = numCells_An5(faces_An5(faceCorresp(:,2),:));
tmp3 = tmp1 - tmp2;
%   ### The face mobility should be associated with delta_G1 - delta_G2, ###
%   ### because a face moves the most when onegrain grow one grain shrink ###
faceMob_dV = abs(tmp3(:,1) - tmp3(:,2));
clear tmp1 tmp2 tmp3

% --- get the Rodrigues Vector corresponding to the Faces ---
[RFvecs_An4, misAs_An4]  = getFaceRFvecs(file_An4, faces_An4);
[RFvecs_An5, misAs_An5]  = getFaceRFvecs(file_An5, faces_An5);



set(0,'defaultAxesLabelFontSize',1.1)
set(0,'defaultAxesFontSize',19)
figure(1)
scatter(face_curvDiff,face_areaDiff,'filled');
xlabel('Integral Face Curvature Difference')
ylabel('Face Area Difference')
figure(2)
scatter(faceCurves_An4(faceCorresp(:,1),2),face_areaDiff,'filled');
xlabel('Integral Face Curvature in An4')
ylabel('Face Area Difference')
figure(3)
scatter(faceCurves_An5(faceCorresp(:,2),2),face_areaDiff,'filled');
xlabel('Integral Face Curvature in An5')
ylabel('Face Area Difference')
figure(4)
scatter(faceCurves_An4(faceCorresp(:,1),2),faceMob_dV,'filled');
xlabel('Integral Face Curvature in An4')
ylabel('The Associated Volume Change')
figure(5)
scatter(face_curvDiff,faceMob_dV,'filled');
xlabel('Integral Face Curvature Difference')
ylabel('The Associated Volume Change')
figure(6)
scatter(misAs_An4(faceCorresp(:,1)),faceMob_dV,'filled');
xlabel('The Misorientation Angle of Grain Face')
ylabel('The Associated Volume Change')

%%
% ### binByMisA(property, misAs, gridSize) ###
% ### data_grid = [leftBoundaryOfCurrentBin, avg(property), countInBin] ###
data_grid = binByMisA(faceMob_dV, misAs_An4(faceCorresp(:,1)), 0.5);
figure(7)
scatter(data_grid(:,1), data_grid(:,2),'filled');
xlabel('Misorientation in An4')
ylabel('The average dV across face')
% % quick check 
% neighbors1 = getNeighList(1736, numNeigh_An4, neighList_An4);
% ID2 = zeros(length(neighbors1),2);
% for i = 1:length(neighbors1)
%     mask = (lookUp(:,1) == neighbors1(i));
%     if sum(mask) > 0
%         ID2(i,:) = lookUp(mask,:);
%     end
% end
% ID2 = sortrows(ID2,2);
% neighbors2 = getNeighList(1, numNeigh_An5, neighList_An5);
