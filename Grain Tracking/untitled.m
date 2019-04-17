file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
load('look_up_table_an4_an5.mat');

% facelabel = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
% tri_curv =  abs(roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)).';
% tri_area = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5).';
% 
% data_raw = [facelabel.'; tri_curv.'; tri_area.'];
% data_face_tmp = calcFaceCurvature(data_raw);

[tracked_uniqueface_an4_uncrop, tracked_uniqueface_an5_uncrop] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_complete_faces');
% face_itg_curv_an5_uncrop = calcFaceItgCurv(file_an5, tracked_uniqueface_an5_uncrop, 'unique_faces');
% 
% obj_face = tracked_uniqueface_an4_uncrop(face_itg_curv_an4(:,2) == 0, :);
% obj_face_idx = (1:length(tracked_uniqueface_an4_uncrop))';
% obj_face_idx = obj_face_idx(face_itg_curv_an4(:,2) == 0);
% obj_face_itg_curv = face_itg_curv_an4(face_itg_curv_an4(:,2) == 0, :);

%%
not_cropped = (ismember(tracked_uniqueface_an5_uncrop, tracked_uniqueface_an5, 'rows'));
face_itg_curv_an5_init = face_itg_curv_an5_uncrop(not_cropped, :);
diff_crop = face_itg_curv_an5 - face_itg_curv_an5_init;
diff_crop = [(1:length(diff_crop))', diff_crop];
diff_crop = sortrows(diff_crop, 2);
large_areadiff_faces = diff_crop(1:10, 1);
diff_crop = sortrows(diff_crop, 3);
large_curvdiff_faces = diff_crop(1:10, 1);



%%
not_cropped = (ismember(tracked_uniqueface_an5, tracked_uniqueface_an5_true, 'rows'));
diff_fc_nnfc_an4 = diff_fc_nnfc_an4(not_cropped);
diff_fc_nnfc_an5 = diff_fc_nnfc_an5(not_cropped);
corner_diff = corner_diff(not_cropped);




%%
idx = 65;
mask_face = ((data_face_tmp(:,1) == obj_face(idx, 1) & data_face_tmp(:,2) == obj_face(idx, 2)) | ...
            (data_face_tmp(:,1) == obj_face(idx, 2) & data_face_tmp(:,2) == obj_face(idx, 1)));
obj_data_face_tmp = data_face_tmp(mask_face,:)

mask_tri = ((facelabel(:,1) == obj_face(idx, 1) & facelabel(:,2) == obj_face(idx, 2)) | ...
            (facelabel(:,1) == obj_face(idx, 2) & facelabel(:,2) == obj_face(idx, 1)));
        
face_tri_curv = tri_curv(mask_tri);


