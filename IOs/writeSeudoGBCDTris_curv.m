%% #################################### Read Data ####################################
folder = 'C:\Users\Richa\Desktop\xiaoting_data\pseudo_GBPD_curv\';
% foler = 'H:\_190730\pseudo_GBPD_curv\';
% % ----------- inner faces -----------
% file_an4 = 'H:\_190730\D3Ds\an0.dream3d';
% file_an5 = 'H:\_190730\D3Ds\an1to0.dream3d';
% load('H:\_190730\matlab_data\look_up_table_an0_an1crop.mat')
% % disp(calcOriginShiftFromTrackedGrains(file_an4, file_an5, look_up_table))
% load('H:\_190730\pseudo_GBPD\GBCDtris_all.mat', 'an0GBCDtris')
% data = an0GBCDtris;
% clear an0GBCDtris
% fname = [folder, '\an0'];
% % ---------------------------------
% file_an4 = 'H:\_190730\D3Ds\an1with2.dream3d';
% file_an5 = 'H:\_190730\D3Ds\an2to1.dream3d';
% load('H:\_190730\matlab_data\look_up_table_an1crop_an2crop.mat')
% % disp(calcOriginShiftFromTrackedGrains(file_an4, file_an5, look_up_table))
% load('H:\_190730\pseudo_GBPD\GBCDtris_all.mat', 'an1with2GBCDtris')
% data = an1with2GBCDtris;
% clear an1with2GBCDtris
% fname = [folder, 'an1with2'];
% % ---------------------------------
% file_an4 = 'H:\_190730\D3Ds\an2to3.dream3d';
% file_an5 = 'H:\_190730\D3Ds\an3with2.dream3d';
% load('H:\_190730\matlab_data\look_up_table_an2crop_an3crop.mat')
% % disp(calcOriginShiftFromTrackedGrains(file_an4, file_an5, look_up_table))
% load('H:\_190730\pseudo_GBPD\GBCDtris_all.mat', 'an2to3GBCDtris')
% data = an2to3GBCDtris;
% clear an2to3GBCDtris
% fname = [folder, 'an2to3'];
% % ---------------------------------
% file_an4 = 'H:\_190730\D3Ds\an3to4.dream3d';
% file_an5 = 'H:\_190730\D3Ds\an4with3.dream3d';
% load('H:\_190730\matlab_data\look_up_table_an3crop_an4crop.mat')
% % disp(calcOriginShiftFromTrackedGrains(file_an4, file_an5, look_up_table))
% load('H:\_190730\pseudo_GBPD\GBCDtris_all.mat', 'an3to4GBCDtris')
% data = an3to4GBCDtris;
% clear an3to4GBCDtris
% fname = [folder, 'an3to4'];
% % ---------------------------------
file_an4 = 'H:\_190730\D3Ds\An4new6_fixedOrigin_smooth.dream3d';
file_an5 = 'H:\_190730\D3Ds\An5new6_cropToAn4.dream3d';
load('H:\_190730\matlab_data\look_up_table_an4_an5crop.mat')
% disp(calcOriginShiftFromTrackedGrains(file_an4, file_an5, look_up_table))
load('H:\_190730\pseudo_GBPD\GBCDtris_all.mat', 'An4new6fixedOriginsmoothGBCDtris')
data = An4new6fixedOriginsmoothGBCDtris;
clear An4new6fixedOriginsmoothGBCDtris
% fname = 'H:\_190730\pseudo_GBPD\An4new6_fixedOrigin_smooth_amortizedDA';
fname = [folder, 'An4new6_fixedOrigin_smooth_includeTLtris'];
% % ---------------------------------

eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;

% ---------------- Get inner faces ----------------
[tracked_uniqueface_an4_inner, tracked_uniqueface_an5_inner] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_inner_faces');
obj_faces_an4 = tracked_uniqueface_an4_inner;
obj_faces_an5 = tracked_uniqueface_an5_inner;


% ---------------- Calc face area and curvature ----------------
face_tmp_an4 = calcFaceItgCurv(file_an4, obj_faces_an4, 'abs', eps_curv, eps_area, eps_min_ang);
face_tmp_an5 = calcFaceItgCurv(file_an5, obj_faces_an5, 'abs', eps_curv, eps_area, eps_min_ang);
face_area_an4 = face_tmp_an4(:,1);
face_itg_curv_an4 = face_tmp_an4(:,2);
face_area_an5 = face_tmp_an5(:,1);
face_itg_curv_an5 = face_tmp_an5(:,2);
face_area_diff = face_area_an5 - face_area_an4;
face_area_diff_ratio = face_area_diff ./ face_area_an4;
face_itg_curv_diff = face_itg_curv_an5 - face_itg_curv_an4;

% ---------------- Filter extreme faces ----------------
mask = (face_area_an5 > 20 & face_area_an4 > 20 & ... 
    face_area_diff_ratio > -0.9 & face_area_diff_ratio < 10);
face_area_an4 = face_area_an4(mask);
face_area_an5 = face_area_an5(mask);
face_area_diff = face_area_diff(mask);
face_itg_curv_an4 = face_itg_curv_an4(mask);
face_itg_curv_an5 = face_itg_curv_an5(mask);
face_itg_curv_diff = face_itg_curv_diff(mask);
obj_faces_an4 = obj_faces_an4(mask, :);
obj_faces_an5 = obj_faces_an5(mask, :);

% ---------------------------------- Amortize curvature change to each triangle ----------------------------------
% face_itg_curv_amortized_diff = face_itg_curv_diff ./ face_area_an4;
face_itg_curv_amortized_diff = face_itg_curv_an5 ./ face_area_an5 - face_itg_curv_an4 ./ face_area_an4;
face_itg_curv_amortized_an4 = face_itg_curv_an4 ./ face_area_an4;
% -----------------------------------------------------------------------------------------------------------

%%% #################################### Data Prepare, Laplacian Smoosth ####################################
% """
% exclude triple line triangles
% """
tri_nodes = 1 + h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')';
node_types = h5read(file_an4,'/DataContainers/TriangleDataContainer/VertexData/NodeType')';
tri_fl = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
tri_fl = sort(tri_fl, 2);
mask = all(tri_fl >= 0, 2);
tri_fl = tri_fl(mask, :);
tri_nodes = tri_nodes(mask, :);

% tri_node_types = node_types(tri_nodes);
% mask_good_tris = (all(tri_node_types == 2, 2) & all(tri_fl > 0, 2));
mask_good_tris = all(tri_fl > 0, 2);

tri_fl = tri_fl(mask_good_tris, :);
data = data(mask_good_tris, :);


% #################################### Assign Resident Face Values to Triangles ####################################
%  ------------------------------------- assign face values to individual triangles -------------------------------------
tri_amortized_itg_curv_diff = zeros(size(tri_fl, 1), 1) * NaN;
tri_amortized_itg_curv_an4 = zeros(size(tri_fl, 1), 1) * NaN;

for i = 1:length(obj_faces_an4)
    mask = (tri_fl(:, 1) == obj_faces_an4(i, 1) & tri_fl(:, 2) == obj_faces_an4(i, 2));
    tri_amortized_itg_curv_diff(mask) = face_itg_curv_amortized_diff(i);
    tri_amortized_itg_curv_an4(mask) = face_itg_curv_amortized_an4(i);
end


mask_trackedface = ~isnan(tri_amortized_itg_curv_diff);
tri_fl = tri_fl(mask_trackedface, :);
tri_amortized_itg_curv_diff = tri_amortized_itg_curv_diff(mask_trackedface);
tri_amortized_itg_curv_an4 = tri_amortized_itg_curv_an4(mask_trackedface);
% fa_absdiff = abs(tri_amortized_a_diff);
% fa_diff_ratio = tri_amortized_a_diff ./ tri_fa;
% fa_absdiff_ratio = fa_absdiff ./ tri_fa;
% idxes = idxes(mask_trackedface);

data = data(mask_trackedface, :);
% curv = curv(mask_trackedface);

%% #################################### Write Seudo GBCDtris file ####################################
% folder = 'C:\Users\Richa\Desktop\xiaoting_data\pseudo_GBPD_curv\';
% fname = [folder, 'allExcept34_inner'];
% load('C:\Users\Richa\Desktop\xiaoting_data\pseudo_GBPD_curv\an0_amortizedDA_GBCDtris_DfMs.mat', 'data')
% data_an0 = data;
% load('C:\Users\Richa\Desktop\xiaoting_data\pseudo_GBPD_curv\an1with2_amortizedDA_GBCDtris_DfMs.mat', 'data')
% data_an1 = data;
% load('C:\Users\Richa\Desktop\xiaoting_data\pseudo_GBPD_curv\an2to3_amortizedDA_GBCDtris_DfMs.mat', 'data')
% data_an2 = data;
% load('H:\_190730\pseudo_GBPD_amorized_curv_change\An4new6_fixedOrigin_smooth_includeTLtris_GBCDtris_DfMs.mat', 'data')
% data_an4 = data;
% data = [data_an0; data_an1; data_an2; data_an4];
% fa_diff = [fa_diff_an0; fa_diff_an1; fa_diff_an2; fa_diff_an3; fa_diff_an4];
% fa_diff_ratio = [fa_diff_ratio_an0; fa_diff_ratio_an1; fa_diff_ratio_an2; ...
%     fa_diff_ratio_an3; fa_diff_ratio_an4];
% fa_absdiff = [fa_absdiff_an0; fa_absdiff_an1; fa_absdiff_an2; fa_absdiff_an3; fa_absdiff_an4];
% fa_absdiff_ratio = [fa_absdiff_ratio_an0; fa_absdiff_ratio_an1; ...
%     fa_absdiff_ratio_an2; fa_absdiff_ratio_an3; fa_absdiff_ratio_an4];
% fname = 'fullTrack_all_states';

face_in_use = '_amortized_curv_change';
data = [data, tri_amortized_itg_curv_diff];
mask = data(:, end) > 0;
csvwrite([fname, face_in_use, '_pos.csv'], data(mask, :));
csvwrite([fname, face_in_use, '_neg.csv'], data(~mask, :));
csvwrite([fname, face_in_use, '.csv'], data);

save([fname, '_GBCDtris_DfMs.mat'], 'data', 'tri_amortized_itg_curv_an4', 'tri_amortized_itg_curv_diff', 'tri_fl', ...
    'face_itg_curv_an4', 'face_itg_curv_diff', 'face_area_an4', 'face_itg_curv_amortized_an4', 'face_itg_curv_amortized_diff', 'obj_faces_an4', 'obj_faces_an5');




%% #################################### Write an4 Distributions ####################################
% file_an4 = 'H:\_190730\D3Ds\An4new6_fixedOrigin_smooth.dream3d';
% load('H:\_190730\pseudo_GBPD\GBCDtris_all.mat', 'An4new6fixedOriginsmoothGBCDtris')
% fname = 'H:\_190730\pseudo_GBPD\An4new6_fixedOrigin_smooth';
file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
load('/Volumes/DATAS/_190730/pseudo_GBPD/GBCDtris_all.mat', 'An4new6fixedOriginsmoothGBCDtris')
fname = '/Volumes/DATAS/_190730/pseudo_GBPD_curv/An4new6_fixedOrigin_smooth';

data = An4new6fixedOriginsmoothGBCDtris;
clear An4new6fixedOriginsmoothGBCDtris
% fname = 'H:\_190730\pseudo_GBPD\An4new6_fixedOrigin_smooth_amortizedDA';

node_types = h5read(file_an4,'/DataContainers/TriangleDataContainer/VertexData/NodeType')';
tri_fl = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
tri_curv =  abs(roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5))';
tri_fl = sort(tri_fl, 2);
mask = all(tri_fl >= 0, 2);
tri_fl = tri_fl(mask, :);
tri_curv = tri_curv(mask, :);

mask_good_tris = all(tri_fl > 0, 2) & tri_curv < 1;

tri_fl = tri_fl(mask_good_tris, :);
data = data(mask_good_tris, :);
tri_curv = tri_curv(mask_good_tris, :);
% 
% csvwrite([fname, '_GBCD.csv'], data);
csvwrite([fname, '_GBHD.csv'], [data, tri_curv]);















