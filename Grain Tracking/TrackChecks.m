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



%% ##### check the calculation of face curvature #####
file = file_An4;
faces = faces_An4;
idx_face = 1846;

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

FCurvs_An4(idx_face,:)



%% ##### check the calculation of face centroids #####
file = file_An5;
faces = faces_An5;
idx_face = 1830;

GA = faces(idx_face, 1);
GB = faces(idx_face, 2);

mask = ((FL(:,1) == GA & FL(:,2) == GB) | (FL(:,1) == GB & FL(:,2) == GA));
faceNodes = triNodes(mask, :);
faceCoords = NCoords(reshape(faceNodes, [], 1),:);

faceCentroid_test = sum(faceCoords)/length(faceCoords)
faceCentroids(idx_face, :)


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







