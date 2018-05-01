function faceCurves = faceCurvatureForTrack(file, faces)
% ##########################################################################
% * Input
%     - faces([N,2]) and faceCorresp([N,1]) is returned by the function 'TrackFace'
% * Output
%     - faceCurvs = [N,2] = [integralArea, integralCurvature]
%         the data was in the same order as faces so the correspodence 
%         returned by 'TrackFace' can be applied directly. 
% ##########################################################################
% ----------------------- load debug data -----------------------

% ---------------------------------------------------------------
    FL = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
    triCurves =  roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5).';
    triAreas = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5).';
    
    data_raw = [FL.'; triCurves.'; triAreas.'];
    data_face_tmp = calcFaceCurvature(data_raw);
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

    % ### note the sign of triangle curvature was decided by winding, can't add the half faces directly ###
    % ### choose the first grain as reference frame to calc face curvature, calc for all faces so there will be a duplicate  ###
    for i = 1:2:length(data_face)
        data_face(i, 4) = data_face(i, 4) + data_face(i+1, 4);
        data_face(i+1, 4) = data_face(i, 4);
        data_face(i, 5) = data_face(i, 4)*data_face(i, 5) - data_face(i+1, 4)*data_face(i+1, 5);
        data_face(i+1, 5) = data_face(i, 5);
    end
    data_face = sortrows(data_face, 1);

    % ### Note that due to the step of reorder labels (13 lines above), the faceLabel won't match faces. ###
    % ### Instead, will match to sort(faces,2). Anyway, the positions are correct so the order can be implicit. ### 
    % ### sum(data_face(:,2:3) == sort(sortrows(faces(faceCorresp(:,1),:),[1,2]),2)); ### 
    faceCurves = data_face(:,4:5);
%     faceCurves = data_face;
end

