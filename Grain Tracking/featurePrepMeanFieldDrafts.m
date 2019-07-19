% """
% Hypothesis: geometry/topology changes in neighbors of grain faces will drag of the face to move.
% Can't assign sign to curvature because the sign of untracked grain faces would be random.
% """  
% ------------------ Get grain face connections ------------------
% nn_faces = getLeftRightFaceConnections(unique_faces, tls, file)
% load('/Volumes/XIAOTING/Ni/working/190425_Hsmooth_energygrad_an5crop.mat', 'nn_faces_an4', 'nn_faces_an5')
% % % % load('/Volumes/XIAOTING/Ni/working/Grain Tracking/data/190221_mig_sign.mat', 'nn_faces_an4', 'nn_faces_an5')


% % ----------------------- load debug data -----------------------
% file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% file_an5 = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth.dream3d';
% load('look_up_table_an4_an5crop.mat')
% 
% % """
% % Threshold values for determining if a triangle is good. Use in calcFaceToGrainCentroidDist.m
% % """
% eps_curv = 1;
% eps_area = 7;
% eps_min_ang = 10;

% ------------------ calculate area for all grain faces ------------------
[faces_an4_all, faces_an5_all, face_corresp] = trackFace(file_an4, file_an5, look_up_table, 0);
face_tmp_an4 = calcFaceItgCurv(file_an4, faces_an4_all, 'as_given', eps_curv, eps_area, eps_min_ang);
face_tmp_an5 = calcFaceItgCurv(file_an5, faces_an5_all, 'as_given', eps_curv, eps_area, eps_min_ang);


% ------------------ calculate area change of grain faces ------------------
% """
% diff = an5 - an4;
% """
face_corresp = sortrows(face_corresp, 1);
tracked_faces_an4 = face_corresp(:,1);
face_corresp_map = containers.Map(face_corresp(:,1),face_corresp(:,2));
area_diff = - face_tmp_an4(:,1);
itg_abscurv_diff = - face_tmp_an4(:,2);
for i = 1:size(faces_an4_all, 1)
    if ismember(i, tracked_faces_an4)
        area_diff(i) = area_diff(i) + face_tmp_an5(face_corresp_map(i), 1);
        itg_abscurv_diff(i) = itg_abscurv_diff(i) + face_tmp_an5(face_corresp_map(i), 2);
    end
end


%%
% ------------------ keep only unique face data in an4 ------------------
faces_an4_all = sort(faces_an4_all, 2);
[~, idx_unique_an4] = unique(faces_an4_all, 'rows');
faces_an4_all = faces_an4_all(idx_unique_an4, :);
area_diff = area_diff(idx_unique_an4);
itg_abscurv_diff = itg_abscurv_diff(idx_unique_an4);



%%
% ------------------ max & avg area change for the neighboring faces ------------------
% nn_faces_an4 = getLeftRightFaceConnections(tracked_uniqueface_an4, tls_an4, file_an4);
% nn_faces_an5 = getLeftRightFaceConnections(tracked_uniqueface_an5, tls_an5, file_an5);

fnnf_area_maxdec = zeros(size(tracked_uniqueface_an4, 1), 1);
fnnf_area_avgdec = zeros(size(tracked_uniqueface_an4, 1), 1);
fnnf_itgcurv_maxdec = zeros(size(tracked_uniqueface_an4, 1), 1);
fnnf_itgcurv_avgdec = zeros(size(tracked_uniqueface_an4, 1), 1);
for i = 1:size(tracked_uniqueface_an4, 1)
    nnf_an4_left = nn_faces_an4{i, 1};
    mask_lnnf_an4 = ismember(faces_an4_all, nnf_an4_left(:, 1:2), 'rows');
    nnf_area_maxdec_left = min(area_diff(mask_lnnf_an4, :));
    nnf_area_avgdec_left = sum(area_diff(mask_lnnf_an4, :))/sum(mask_lnnf_an4);
    nnf_itgcurv_maxdec_left = min(itg_abscurv_diff(mask_lnnf_an4, :));
    nnf_itgcurv_avgdec_left = sum(itg_abscurv_diff(mask_lnnf_an4, :))/sum(mask_lnnf_an4);
    nnf_an4_right = nn_faces_an4{i, 2};
    mask_rnnf_an4 = ismember(faces_an4_all, nnf_an4_right(:, 1:2), 'rows');
    nnf_area_maxdec_right = min(area_diff(mask_rnnf_an4, :));
    nnf_area_avgdec_right = sum(area_diff(mask_rnnf_an4, :))/sum(mask_rnnf_an4);
    nnf_itgcurv_maxdec_right = min(itg_abscurv_diff(mask_rnnf_an4, :));
    nnf_itgcurv_avgdec_right = sum(itg_abscurv_diff(mask_rnnf_an4, :))/sum(mask_rnnf_an4);
    
    fnnf_area_maxdec(i) = nnf_area_maxdec_right - nnf_area_maxdec_left;
    fnnf_area_avgdec(i) = nnf_area_avgdec_right - nnf_area_avgdec_left;
    fnnf_itgcurv_maxdec(i) = nnf_itgcurv_maxdec_right - nnf_itgcurv_maxdec_left;
    fnnf_itgcurv_avgdec(i) = nnf_itgcurv_avgdec_right - nnf_itgcurv_avgdec_left;
end


%%
% ------------------ fraction of neighboring faces: twin & positive curvature ------------------
% ##### clean data_face to be for complete grain faces #####
mask_reverse = (data_face_an4(:, 1) > data_face_an4(:, 2));
data_face_an4(mask_reverse, 4) = - data_face_an4(mask_reverse, 4);
data_face_an4(:, 1:2) = sort(data_face_an4(:, 1:2), 2);
faces_unique_an4 = sortrows(data_face_an4, [1,2]);
cnt = 1;
while cnt < length(faces_unique_an4)
    if faces_unique_an4(cnt, 1) == faces_unique_an4(cnt+1, 1) && faces_unique_an4(cnt, 2) == faces_unique_an4(cnt+1, 2)
        face_area = faces_unique_an4(cnt, 3) + faces_unique_an4(cnt+1, 3);
        faces_unique_an4(cnt, 4) = (faces_unique_an4(cnt, 3)*faces_unique_an4(cnt, 4) ...
                            + faces_unique_an4(cnt+1, 3)*faces_unique_an4(cnt+1, 4)) / face_area;
        faces_unique_an4(cnt, 3) = face_area;
        faces_unique_an4(cnt+1, :) = [];
    end
    cnt = cnt + 1;
end 

% ##### identify twin label for faces_unique #####
tri_istwin_an4 = boolean(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/TwinBoundary'))';
fl_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
mask_an4 = all(fl_an4>0, 2);
fl_an4 = fl_an4(mask_an4, :);
tri_istwin_an4 = tri_istwin_an4(mask_an4, :);
fl_an4 = sort(fl_an4, 2);

face_unique_twins_an4 = zeros(length(faces_unique_an4), 1);
for i = 1:length(faces_unique_an4)
    mask_an4 = (fl_an4(:,1) == faces_unique_an4(i, 1) & fl_an4(:,2) == faces_unique_an4(i, 2));
    if all(tri_istwin_an4(mask_an4))
        face_unique_twins_an4(i) = true;
    else
        if any(tri_istwin_an4(mask_an4))
            warning(['D3D problem, an4, pair ', num2str(i)]);
        end
    end
end

% ##### calc fractions #####
neighfrac_pos_an4 = zeros(length(nn_faces_an4), 1);
twinfrac_pos_an4 = zeros(length(nn_faces_an4), 1);

for i = 1:size(tracked_uniqueface_an4, 1)
    % ----- left connections -----
    nnf_an4_left = nn_faces_an4{i, 1};
    mask_lnnf_an4 = ismember(faces_unique_an4(:, 1:2), nnf_an4_left(:, 1:2), 'rows');
    twins_lnnf = face_unique_twins_an4(mask_lnnf_an4);
    tmp = faces_unique_an4(mask_lnnf_an4, :);
    lnnf_curv_an4 = tmp(:,3);
    mask_reverse = (tmp(:,1) == tracked_uniqueface_an4(i, 1));
    lnnf_curv_an4(mask_reverse) = -lnnf_curv_an4(mask_reverse);
    
    % ----- right connections -----
    nnf_an4_right = nn_faces_an4{i, 2};
    mask_rnnf_an4 = ismember(faces_unique_an4(:, 1:2), nnf_an4_right(:, 1:2), 'rows');
    twins_rnnf = face_unique_twins_an4(mask_rnnf_an4);
    tmp = faces_unique_an4(mask_rnnf_an4, :);
    rnnf_curv_an4 = tmp(:,3);
    mask_reverse = (tmp(:,1) == tracked_uniqueface_an4(i, 2));
    rnnf_curv_an4(mask_reverse) = -rnnf_curv_an4(mask_reverse);
    
    % ----- calc fractions -----
    num_neighs = (length(lnnf_curv_an4) + length(rnnf_curv_an4));
    neighfrac_pos_an4(i) = (sum(lnnf_curv_an4 > 0) + sum(rnnf_curv_an4 > 0 )) / num_neighs;
    twinfrac_pos_an4(i) = (sum(twins_rnnf) + sum(twins_lnnf)) / num_neighs;
end




