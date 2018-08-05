data = [faceCorresp, FA_diff, FCurv_diff];
data = sortrows(data, [4, 3]);

faceMap_An4 = containers.Map(faceCorresp(:,1), num2cell(faces_An4(faceCorresp(:,1),:),2));
faceMap_An5 = containers.Map(faceCorresp(:,2), num2cell(faces_An5(faceCorresp(:,2),:),2));


%% ##### The Complete Tracked Grains #####
SG_An4 = h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
SG_An5 = h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
SG_An4(1) = [];
SG_An5(1) = [];

mask_lookUp = logical([SG_An4(lookUp(:,1)), SG_An5(lookUp(:,2))]);
lookUp_CG = lookUp(all(~mask_lookUp, 2), :);


%% ##### The Grains that growed from being inner to touch surface #####
[faces_An4, faces_An5, faceCorresp_AllF] = TrackFace(file_An4, file_An5, lookUp, false);
faceInfo_all = [faces_An4(faceCorresp_AllF(:,1),:), faces_An5(faceCorresp_AllF(:,2),:)];

SG_An4 = h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
SG_An5 = h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
SG_An4(1) = [];
SG_An5(1) = [];
SGinfo_all = [SG_An4(faceInfo_all(:,1:2)), SG_An5(faceInfo_all(:,3:4))];

% tmp1 = (SGinfo_all(:,1)==1 & SGinfo_all(:,2)==1 & SGinfo_all(:,3)==0 & SGinfo_all(:,4)==0);
tmp2 = (SGinfo_all(:,1)==0 & SGinfo_all(:,2)==0 & SGinfo_all(:,3)==1 & SGinfo_all(:,4)==1);
mask_growToSurf = tmp2;

% FItgCurvs_An4 = faceCurvatureForTrack(file_An4, faces_An4);
% FItgCurvs_An5 = faceCurvatureForTrack(file_An5, faces_An5);
diff_tmp = FItgCurvs_An5(faceCorresp_AllF(:,2),:) - FItgCurvs_An4(faceCorresp_AllF(:,1),:);
FA_diff = diff_tmp(mask_growToSurf,1);


%% ##### Plot Triangle Curvature Distribution #####
set(0,'defaultAxesLabelFontSize',1.1)
set(0,'defaultAxesFontSize',19)

FL_An4 = double(h5read(file_An4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
triCurves_An4 =  abs(roundn(h5read(file_An4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)).';
mask_An4 = ~(any(FL_An4 <= 0, 2) | triCurves_An4 > 100);
triCurves_An4 = triCurves_An4(mask_An4);
FL_An5 = double(h5read(file_An5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
triCurves_An5 =  abs(roundn(h5read(file_An5,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)).';
mask_An5 = ~(any(FL_An5 <= 0, 2) | triCurves_An5 > 100);
triCurves_An5 = triCurves_An5(mask_An5);

histogram(triCurves_An4,100)
hold on
histogram(triCurves_An5,100)
set(gca, 'YScale', 'log')
aveTriCurv_An4 = sum(triCurves_An4)/length(triCurves_An4);
aveTriCurv_An5 = sum(triCurves_An5)/length(triCurves_An4);
legend(['aveTriCurv, An4 = ', num2str(aveTriCurv_An4)], ['aveTriCurv, An5 = ', num2str(aveTriCurv_An5)])
xlabel('Triangle Curvature, \mum^{-1}')
ylabel('# Triangles')


%% ##### Plot Triangle Curvature Distribution #####
% get the coords of one face pair

idx = 1;
idx_An4 = faceCorresp(idx, 1);
idx_An5 = faceCorresp(idx, 2);
label_An4 = faces_An4(idx_An4, :);
label_An5 = faces_An5(idx_An5, :);

FL_An4 = double(h5read(file_An4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
triNodes_An4 = 1 + double(h5read(file_An4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
NCoords_An4 = double(h5read(file_An4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
FL_An5 = double(h5read(file_An5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
triNodes_An5 = 1 + double(h5read(file_An5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
NCoords_An5 = double(h5read(file_An5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

mask_An4 = ((FL_An4(:,1) == label_An4(1) & FL_An4(:,2) == label_An4(2)) | (FL_An4(:,1) == label_An4(2) & FL_An4(:,2) == label_An4(1)));
mask_An5 = ((FL_An5(:,1) == label_An5(1) & FL_An5(:,2) == label_An5(2)) | (FL_An5(:,1) == label_An5(2) & FL_An5(:,2) == label_An5(1)));

thisFaceNodes_An4 = triNodes_An4(mask_An4,:);
thisFaceNodes_An5 = triNodes_An5(mask_An5,:);
thisFaceNodes_An4 = reshape(thisFaceNodes_An4, [], 1);
thisFaceNodes_An5 = reshape(thisFaceNodes_An5, [], 1);
triNodeCoords_An4 = NCoords_An4(thisFaceNodes_An4,:);
triNodeCoords_An5 = NCoords_An5(thisFaceNodes_An5,:);

tmpNodes_An4 = triNodes_An4(mask_An4,:);
tmpNodes_An5 = triNodes_An5(mask_An5,:);
%%
trisurf( tmpNodes_An4, NCoords_An4(:,1), NCoords_An4(:,2), NCoords_An4(:,3) );
rotate3d on
hold on
trisurf( tmpNodes_An5, NCoords_An5(:,1), NCoords_An5(:,2), NCoords_An5(:,3),'Facecolor',[218/255 168/255 32/255] );
rotate3d on




%% ##### Find the faceID of a face between two grains #####
[~, ~, faceCorresp_CF] = TrackFace(file_An4, file_An5, lookUp, true);

faceInfo_CF = [C_An4, faces_An4(faceCorresp_CF(:,1),:)];
faceInfo_CF = sortrows(faceInfo_CF);

faceIDs = h5read(file_An4, '/DataContainers/TriangleDataContainer/FaceData/FeatureFaceId')';
FL = h5read(file_An4, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';

GA = faceInfo_CF(68,2)
GB = faceInfo_CF(68,3)
mask = ((FL(:,1) == GA & FL(:,2) == GB) | (FL(:,1) == GB & FL(:,2) == GA));
unique(faceIDs(mask))




clc


file_An4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_mesh.dream3d');
file_An5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh.dream3d');

FL_An4 = h5read(file_An4, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
FL_An5 = h5read(file_An5, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
faceId_An4 = h5read(file_An4, '/DataContainers/TriangleDataContainer/FaceData/FeatureFaceId')';
faceId_An5 = h5read(file_An5, '/DataContainers/TriangleDataContainer/FaceData/FeatureFaceId')';


tmp = [(1:length(FCentrs_diff))', FCentrs_diff];
sortTmp = sortrows(tmp, -2);


objFace = 1;
objLabel_An5 = trackedFL_An5(sortTmp(objFace, 1), :);
mask = (FL_An5(:,1) == objLabel_An5(1) & FL_An5(:,2) == objLabel_An5(2));
unique(FL_An5(mask, :), 'rows')
unique(faceId_An5(mask))

objLabel_An4 = trackedFL_An4(sortTmp(objFace, 1), :);
mask = (FL_An4(:,1) == objLabel_An4(1) & FL_An4(:,2) == objLabel_An4(2));
unique(FL_An4(mask, :), 'rows')
unique(faceId_An4(mask))



%% ##### get FeatureFaceId for the faces whose centroids moved the most #####
clc

file_An4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_mesh.dream3d');
file_An5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh.dream3d');

faceId_An4 = h5read(file_An4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
faceId_An5 = h5read(file_An5, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
nTris_An4 = h5read(file_An4, '/DataContainers/TriangleDataContainer/FaceFeatureData/NumTriangles')';
nTris_An5 = h5read(file_An5, '/DataContainers/TriangleDataContainer/FaceFeatureData/NumTriangles')';
SG_An4 = h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
SG_An5 = h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
faceId_An4(1,:) = [];
faceId_An5(1,:) = [];
faceId_An4 = sort(faceId_An4, 2);
faceId_An5 = sort(faceId_An5, 2);
SG_An4(1) = [];
SG_An5(1) = [];

tmp = [(1:length(FCentrs_diff))', FCentrs_diff];
sortTmp = sortrows(tmp, -2);

cond = true;
while cond
    objFace = round(rand()*14000);
    trackedFL_An4 = sort(trackedFL_An4, 2);

    objLabel_An4 = trackedFL_An4(sortTmp(objFace, 1), :);
    [idx_An4, ~] = find(faceId_An4(:,1) == objLabel_An4(1) & faceId_An4(:,2) == objLabel_An4(2));
    objLabel_An5 = trackedFL_An5(sortTmp(objFace, 1), :);
    [idx_An5, ~] = find((faceId_An5(:,1) == objLabel_An5(1) & faceId_An5(:,2) == objLabel_An5(2)) | (faceId_An5(:,1) == objLabel_An5(2) & faceId_An5(:,2) == objLabel_An5(1)));
    if ~isempty(faceId_An4(idx_An4, :)) && ~isempty(faceId_An5(idx_An5, :))
        cond1 = any([SG_An4(faceId_An4(idx_An4, :)); SG_An5(faceId_An5(idx_An5, :))]==1);
        cond2 = (nTris_An4(idx_An4)<300 || nTris_An5(idx_An5)<300);
        cond = (cond1 || cond2);
    end
end
disp(['FaceId in An4:   ', num2str(idx_An4), ',   size=', num2str(nTris_An4(idx_An4))])
disp(['FaceLabel in An4:  [', num2str(faceId_An4(idx_An4, :)), ']', ]);
disp(['FaceId in An5:   ', num2str(idx_An5), ',   size=', num2str(nTris_An5(idx_An5))])
disp(['FaceLabel in An5:  [', num2str(faceId_An5(idx_An5, :)), ']', ]);



