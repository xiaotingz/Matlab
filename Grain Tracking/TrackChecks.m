%% ##### check the track of faces #####
faces_An4_tracked = faces_An4(faceCorresp(:,1), :);
faces_An5_tracked = faces_An5(faceCorresp(:,2), :);
map = containers.Map(lookUp(:,1), lookUp(:,2));
faces_An5_trackedFromAn4 = zeros(size(faceCorresp));
for i = 1:length(faceCorresp)
    faces_An5_trackedFromAn4(i,1) = map(faces_An4_tracked(i,1));
    faces_An5_trackedFromAn4(i,2) = map(faces_An4_tracked(i,2));
end
sum(sum(faces_An5_tracked ~= faces_An5_trackedFromAn4))



%% ##### check if really only the complete grain faces are used #####
[faces_An4, faces_An5, faceCorresp_AllF] = TrackFace(file_An4, file_An5, lookUp, false);
[~, ~, faceCorresp_CF] = TrackFace(file_An4, file_An5, lookUp, true);

faceInfo = [faces_An4(faceCorresp_AllF(:,1),:), faces_An5(faceCorresp_AllF(:,2),:)];
faceInfo_CF = [faces_An4(faceCorresp_CF(:,1),:), faces_An5(faceCorresp_CF(:,2),:)];

SG_An4 = h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
SG_An5 = h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
SG_An4(1) = [];
SG_An5(1) = [];
SGinfo_all = [SG_An4(faceInfo(:,1:2)), SG_An5(faceInfo(:,3:4))];
SGinfo_CF = [SG_An4(faceInfo_CF(:,1:2)), SG_An5(faceInfo_CF(:,3:4))];

tmp = [any(SGinfo_all(:,1:2)==0,2), any(SGinfo_all(:,3:4)==0,2)];
if sum(all(tmp==1, 2)) ~= length(faceCorresp_CF)
    warning('HEY, completeness of the matched faces are not right!')
end

tmp1 = (SGinfo_all(:,1)==1 & SGinfo_all(:,2)==1 & SGinfo_all(:,3)==0 & SGinfo_all(:,4)==0);
tmp2 = (SGinfo_all(:,1)==0 & SGinfo_all(:,2)==0 & SGinfo_all(:,3)==1 & SGinfo_all(:,4)==1);
mask_growToSurf = (tmp1 | tmp2);


%% ##### check the calculation of face curvature #####
file = file_An4;
faces = faces_An4;
FItgCurvs = FItgCurvs_An4;
idx_face = 280;

FL = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
triCurves =  abs(roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)).';
triAreas = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5).';
data = [FL, triCurves, triAreas];

GA = faces(idx_face, 1);
GB = faces(idx_face, 2);
mask_a = (data(:,1) == GA & data(:,2) == GB);
mask_b = (data(:,1) == GB & data(:,2) == GA);
data_select_a = data(mask_a,:);
data_select_b = data(mask_b,:);

faceA = sum(data_select_a(:,4)) + sum(data_select_b(:,4))
faceIC = sum(data_select_a(:,3).*data_select_a(:,4)) + sum(data_select_b(:,3).*data_select_b(:,4))

FItgCurvs(idx_face,:)



%% ##### check the calculation of face centroids #####
file = file_An4;
faces = faces_An4;
FCentrs = FCentrs_An4;

FL = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
triNodes = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
NCoords = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
idx_face = 23171;

GA = faces(idx_face, 1);
GB = faces(idx_face, 2);

mask = ((FL(:,1) == GA & FL(:,2) == GB) | (FL(:,1) == GB & FL(:,2) == GA));
faceNodes = triNodes(mask, :);
faceCoords = NCoords(reshape(faceNodes, [], 1),:);

faceCentroid_test = sum(faceCoords)/length(faceCoords)
FCentrs(idx_face, :)


%% ##### check TrackUniqueFace #####
UFcorresp = TrackUniqueFace(UFlabels_An4, UFlabels_An5, lookUp);
GIDMap_4to5 = containers.Map(lookUp(:,1), num2cell(lookUp(:,2)));

% --- check the last element and some random elements ---
for i = 1:5
    if i == 1
        idx = length(UFcorresp);
    else
        idx = randi(length(UFcorresp));
    end
    idx_An4 = UFcorresp(idx,1);
    idx_An5 = UFcorresp(idx,2);
    label_An4 = UFlabels_An4(idx_An4, :);
    label_An4in5 = [GIDMap_4to5(label_An4(1)), GIDMap_4to5(label_An4(2))]
    label_An5 = UFlabels_An5(idx_An5, :)
    if ~ ((label_An4in5(1) == label_An5(1) && label_An4in5(2) == label_An5(2)) || (label_An4in5(1) == label_An5(2) && label_An4in5(2) == label_An5(1)))
        warning(['NOT matching at face=', num2str(idx)])
    end
end







