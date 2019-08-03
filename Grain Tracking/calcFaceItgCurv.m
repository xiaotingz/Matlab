function face_itgcurv = calcFaceItgCurv(file, faces, sign_req, eps_curv, eps_area, eps_min_ang)
% ##########################################################################
% * Input
%     - faces([n,2]) and faceCorresp([n,1]) 
%         is returned by the function trackFace.m or trackUniqueFace.m
%     - sign_req = 'signed_resident_left' | 'abs'
%         'signed_resident_left': calculates signed face_itg_curv assuming resident grain on left.
%         'abs': calculate face_itg_abs_curv
%     - eps_
%         thresholds for good quality triangles. Especially important for data by Hsmooth.
% * Output
%     - faceCurvs = [n,2] = [integralArea, integralAbsCurvature]
%         the data was in the same order as faces so the correspodence 
%         returned by 'TrackFace' can be applied directly. 
% * NOTE
%     - face_itgcurv contain duplicates of entire face data, not half face data. 
%     - integralCurvature = sum(triArea * abs(triCurvature)), direction can't be told without the grain of interest.
%     - Dependency: calcFaceCurvature.m
% ##########################################################################
% ----------------------- load debug data -----------------------
% file = file_an4;
% faces = faces_an4;
% ---------------------------------------------------------------
facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
tri_curv =  roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
tri_area = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
tri_min_ang = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';

if strcmp(sign_req, 'abs')
    tri_curv = abs(tri_curv);
else
    if ~ strcmp(sign_req, 'signed_resident_left')
        warning('Wrong input sign_req, use ''signed_resident_left'' or ''abs''')
        return
    end
end
    
% data_raw = [facelabel.'; tri_curv.'; tri_area.'];
data_raw = [facelabel'; tri_curv'; tri_area'; tri_min_ang'];

% """
% data_face_tmp = [label_1, label_2, area, itg_curv/area]
% """
data_face_tmp = calcFaceCurvature(data_raw, eps_curv, eps_area, eps_min_ang);

% if strcmp(input, 'all_faces')
%     --------------------------------------------------------------------------------------
%     """ BUG INSIDE """
%     % ### then make data_face the same format as faces, so that the indexes can be used ###
%     
%     data_face = zeros(length(faces),4);
%     data_face(:,1:2) = faces;
%     idx = 1;
%     for i = 1 : length(data_face)
%         if (data_face(i,1) == data_face_tmp(idx,1)) && (data_face(i,2) == data_face_tmp(idx,2))
%             data_face(i,:) = data_face_tmp(idx, :);
%             idx = idx + 1;
%         end
%     end
%     data_face = [(1:length(data_face))', data_face]; 
% 
%     % ### reorder the face labels to make find unique faces easier ###
%     data_face(:,2:3) = sort(data_face(:,2:3),2);
%     data_face = sortrows(data_face, [2,3]);
%     for i = 1:2:length(data_face)
%         data_face(i, 4) = data_face(i, 4) + data_face(i+1, 4);
%         data_face(i, 5) = data_face(i, 4)*data_face(i, 5) + data_face(i+1, 4)*data_face(i+1, 5);
%         data_face(i+1, 4) = data_face(i, 4);
%         data_face(i+1, 5) = data_face(i, 5);
%     end
% 
%     % ### sort the order back ###
%     data_face = sortrows(data_face, 1);
% 
%     % ### Note that due to the step of reorder labels (13 lines above), the faceLabel won't match faces. ###
%     % ### Instead, will match to sort(faces,2). Anyway, the positions are correct so the order can be implicit. ### 
%     % ### sum(data_face(:,2:3) == sort(sortrows(faces(faceCorresp(:,1),:),[1,2]),2)); ### 
%     face_itgcurv = data_face(:,4:5);
%     %     faceCurves = data_face;
%     --------------------------------------------------------------------------------------

face_itgcurv = zeros(size(faces));
for i = 1:size(faces, 1)
    mask_1 = data_face_tmp(:,1) == faces(i,1) & data_face_tmp(:,2) == faces(i,2);
    mask_2 = data_face_tmp(:,1) == faces(i,2) & data_face_tmp(:,2) == faces(i,1);
    mask_face = (mask_1 | mask_2);
    
    face_itgcurv(i, 1) = sum(data_face_tmp(mask_face, 3));
    
    if strcmp(sign_req, 'signed_resident_left')
        face_itgcurv(i, 2) = - sum(data_face_tmp(mask_1, 3) .* data_face_tmp(mask_1, 4)) ...
                            + sum(data_face_tmp(mask_2, 3) .* data_face_tmp(mask_2, 4));
    elseif strcmp(sign_req, 'abs')
        face_itgcurv(i, 2) = sum(data_face_tmp(mask_face, 3) .* data_face_tmp(mask_face, 4));
    end
    
end

end



% %%
% % ############################# CHECK #############################
% file = file_an4;
% faces = tracked_uniqueface_an4;
% face_area = face_area_an4;
% face_avg_cur = face_avg_curv_an4;
% face_itg_curv = face_itg_curv_an4;
%  
% 
% facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
% tri_curv =  roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
% tri_area = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
% tri_min_ang = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
% mask = all(facelabel>0, 2) & abs(tri_curv)<eps_acurv & tri_area<eps_area & tri_min_ang>eps_min_ang;
% facelabel = facelabel(mask, :);
% tri_curv = tri_curv(mask);
% tri_area = tri_area(mask);
% tri_min_ang = tri_min_ang(mask);
% %%
% rng('shuffle')
% idx = randi(length(faces));
% 
% mask_1 = (facelabel(:,1)==faces(idx, 1) & facelabel(:,2)==faces(idx, 2));
% mask_2 = (facelabel(:,1)==faces(idx, 2) & facelabel(:,2)==faces(idx, 1));
% mask = mask_1 | mask_2;
% itg_abscurv = sum(tri_area(mask) .* abs(tri_curv(mask)));
% itg_curv_left = - sum(tri_area(mask_1) .* tri_curv(mask_1)) + sum(tri_area(mask_2) .* tri_curv(mask_2));
% avg_curv = sum(tri_curv(mask)) / sum(mask);
% avg_abscurv = sum(abs(tri_curv(mask))) / sum(mask);
% area = sum(tri_area(mask));
% 
% disp(['pair = ', num2str(idx)]);
% disp(['calced_area = ', num2str(area), ';   check_area = ', num2str(face_area(idx))]);
% disp(['calced_avg_curv = ', num2str(avg_curv)]);
% disp(['calced_avg_abscurv = ', num2str(avg_abscurv), ';   check_avg_abscurv = ', num2str(face_avg_cur(idx))]);
% disp(['calced_itg_abscurv = ', num2str(itg_abscurv), ';   check_itg_abscurv = ', num2str(face_itg_curv(idx))]);
% % disp(['calced_itg_curv_left = ', num2str(itg_curv_left), ';   check_itg_curv_left = ', num2str(face_itgcurv_an4_left(idx))]);
% disp(' ')


























