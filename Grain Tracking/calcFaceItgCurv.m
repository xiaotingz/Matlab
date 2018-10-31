function face_itgcurv = calcFaceItgCurv(file, faces, input)
% ##########################################################################
% * Input
%     - faces([N,2]) and faceCorresp([N,1]) is returned by the function
%     trackFace.m or trackUniqueFace.m
% * Output
%     - faceCurvs = [N,2] = [integralArea, integralCurvature]
%         the data was in the same order as faces so the correspodence 
%         returned by 'TrackFace' can be applied directly. 
% * NOTE
%     - integralCurvature = sum(triArea * abs(triCurvature)), direction can't be told without the grain of interest.
% ##########################################################################
% ----------------------- load debug data -----------------------
% file = file_an4;
% faces = faces_an4;
% ---------------------------------------------------------------
FL = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
triCurves =  abs(roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)).';
triAreas = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5).';

data_raw = [FL.'; triCurves.'; triAreas.'];
% """
% data_face_tmp = [label_1, label_2, area, itg_curv/area]
% """
data_face_tmp = calcFaceCurvature(data_raw);

if strcmp(input, 'all_faces')
    % ### then make data_face the same format as faces, so that the indexes can be used ###
    data_face = zeros(length(faces),4);
    data_face(:,1:2) = faces;
    idx = 1;
    for i = 1 : length(data_face)
        if (data_face(i,1) == data_face_tmp(idx,1)) && (data_face(i,2) == data_face_tmp(idx,2))
            data_face(i,:) = data_face_tmp(idx, :);
            idx = idx + 1;
        end
    end
    data_face = [(1:length(data_face))', data_face]; 

    % ### reorder the face labels to make find unique faces easier ###
    data_face(:,2:3) = sort(data_face(:,2:3),2);
    data_face = sortrows(data_face, [2,3]);
    for i = 1:2:length(data_face)
        data_face(i, 4) = data_face(i, 4) + data_face(i+1, 4);
        data_face(i, 5) = data_face(i, 4)*data_face(i, 5) + data_face(i+1, 4)*data_face(i+1, 5);
        data_face(i+1, 4) = data_face(i, 4);
        data_face(i+1, 5) = data_face(i, 5);
    end

    % ### sort the order back ###
    data_face = sortrows(data_face, 1);

    % ### Note that due to the step of reorder labels (13 lines above), the faceLabel won't match faces. ###
    % ### Instead, will match to sort(faces,2). Anyway, the positions are correct so the order can be implicit. ### 
    % ### sum(data_face(:,2:3) == sort(sortrows(faces(faceCorresp(:,1),:),[1,2]),2)); ### 
    face_itgcurv = data_face(:,4:5);
    %     faceCurves = data_face;
elseif strcmp(input, 'unique_faces')
    face_itgcurv = zeros(size(faces));
    for i = 1:length(faces)
        mask_face = ((data_face_tmp(:,1) == faces(i,1) & data_face_tmp(:,2) == faces(i,2)) | ...
            (data_face_tmp(:,1) == faces(i,2) & data_face_tmp(:,2) == faces(i,1)));
        face_itgcurv(i, 1) = sum(data_face_tmp(mask_face, 3));
        face_itgcurv(i, 2) = sum(data_face_tmp(mask_face, 3).*data_face_tmp(mask_face, 4));
    end
else
    warning('The third argument of calcFaceItgCurv should be either all_faces or unique_faces.')
end

end
