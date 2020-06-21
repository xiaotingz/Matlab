% This file finds the grain faces that are truely connected to the sample
% free surface. In the context of HEDM pillars, those are the faces with
% facelabel=0

%% #################################### Read Data ####################################
% % ----------- inner faces -----------
% file_an4 = 'H:\_190730\D3Ds\an0.dream3d';
% file_an5 = 'H:\_190730\D3Ds\an1to0.dream3d';
% load('H:\_190730\matlab_data\look_up_table_an0_an1crop.mat')
% fname = 'H:\_190730\pseudo_GBPD\an0_amortizedDA';
% % ---------------------------------
% file_an4 = 'H:\_190730\D3Ds\an1with2.dream3d';
% file_an5 = 'H:\_190730\D3Ds\an2to1.dream3d';
% load('H:\_190730\matlab_data\look_up_table_an1crop_an2crop.mat')
% fname = 'H:\_190730\pseudo_GBPD\an1with2_amortizedDA';
% % ---------------------------------
% file_an4 = 'H:\_190730\D3Ds\an2to3.dream3d';
% file_an5 = 'H:\_190730\D3Ds\an3with2.dream3d';
% load('H:\_190730\matlab_data\look_up_table_an2crop_an3crop.mat')
% fname = 'H:\_190730\pseudo_GBPD\an2to3_amortizedDA';
% % % % ---------------------------------
% % % file_an4 = 'H:\_190730\D3Ds\an3to4.dream3d';
% % % file_an5 = 'H:\_190730\D3Ds\an4with3.dream3d';
% % % load('H:\_190730\matlab_data\look_up_table_an3crop_an4crop.mat')
% % % fname = 'H:\_190730\pseudo_GBPD\an3to4_amortizedDA';
% % ---------------------------------
file_an4 = 'H:\_190730\D3Ds\An4new6_fixedOrigin_smooth.dream3d';
file_an5 = 'H:\_190730\D3Ds\An5new6_cropToAn4.dream3d';
load('H:\_190730\matlab_data\look_up_table_an4_an5crop.mat')
fname = 'H:\_190730\pseudo_GBPD\An4new6_fixedOrigin_smooth_amortizedDA';
% % ---------------------------------

load([fname, '_GBCDtris.mat'], 'data', 'tri_fa', 'tri_amortized_a_diff', 'tri_fl', ...
     'face_area_amortized_diff', 'obj_faces_an4', 'obj_faces_an5');
[tracked_uniqueface_an4_inner, tracked_uniqueface_an5_inner] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_inner_faces');
mask_good_faces = ismember(tracked_uniqueface_an4_inner, obj_faces_an4, 'rows');
obj_faces_an5 = obj_faces_an5(mask_good_faces, :);

% -------- Idenfity true SF grains --------
face_featureid_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels').';
face_featureid_an4(1,:) = [];
mask_true_surf_g_an4 = any(face_featureid_an4==0, 2);
true_surf_g_idx_an4 = unique(face_featureid_an4(mask_true_surf_g_an4, :));
true_surf_g_idx_an4(true_surf_g_idx_an4<=0) = [];
face_featureid_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels').';
face_featureid_an5(1,:) = [];
mask_true_surf_g_an5 = any(face_featureid_an5==0, 2);
true_surf_g_idx_an5 = unique(face_featureid_an5(mask_true_surf_g_an5, :));
true_surf_g_idx_an5(true_surf_g_idx_an5<=0) = [];

mask_touch_fs_f = sum(ismember(obj_faces_an4, true_surf_g_idx_an4), 2) == 2 | ...
    sum(ismember(obj_faces_an5, true_surf_g_idx_an5), 2) == 2;
faces_touch_fs = obj_faces_an4(mask_touch_fs_f, :);
mask_touch_fs_tri = ismember(tri_fl, faces_touch_fs, 'rows');

face_in_use = '_touchSF';
data_fs = data(mask_touch_fs_tri, :);
mask = data_fs(:, end) > 0;
csvwrite([fname, face_in_use, '_pos.csv'],data_fs(mask, :));
csvwrite([fname, face_in_use, '_neg.csv'],data_fs(~mask, :));
csvwrite([fname, face_in_use, '.csv'],data_fs);

face_in_use = '_noSF';
data_no_fs = data(~mask_touch_fs_tri, :);
mask = data_no_fs(:, end) > 0;
csvwrite([fname, face_in_use, '_pos.csv'],data_no_fs(mask, :));
csvwrite([fname, face_in_use, '_neg.csv'],data_no_fs(~mask, :));
csvwrite([fname, face_in_use, '.csv'],data_no_fs);

save([fname, '_maskTouchFS.mat'], 'mask_touch_fs_tri')

% file = 'GBCDtris_noSmallFace_A.txt';
% fileID = fopen(file,'w');
% format = '%7.4f    %7.4f    %7.4f    %7.4f    %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f   %7.4f\n';
% for i = 1:length(data)
%     fprintf(fileID,format,data(i, 1),data(i, 2),data(i, 3),data(i, 4),data(i, 5),data(i, 6), ...
%         data(i, 7),data(i, 8),data(i, 9), data(i, 10), fa(i));
% end
% fclose('all');



%% #################################### Concat all states ####################################
clear

load('H:\_190730\pseudo_GBPD\an0_amortizedDA_GBCDtris.mat', 'data')
load('H:\_190730\pseudo_GBPD\an0_amortizedDA_maskTouchFS.mat')
data_all = data;
mask_trouch_fs_tri_all = mask_touch_fs_tri;

load('H:\_190730\pseudo_GBPD\an1with2_amortizedDA_GBCDtris.mat', 'data')
load('H:\_190730\pseudo_GBPD\an1with2_amortizedDA_maskTouchFS.mat')
data_all = [data_all; data];
mask_trouch_fs_tri_all = [mask_trouch_fs_tri_all; mask_touch_fs_tri];

load('H:\_190730\pseudo_GBPD\an2to3_amortizedDA_GBCDtris.mat', 'data')
load('H:\_190730\pseudo_GBPD\an2to3_amortizedDA_maskTouchFS.mat')
data_all = [data_all; data];
mask_trouch_fs_tri_all = [mask_trouch_fs_tri_all; mask_touch_fs_tri];

% load('H:\_190730\pseudo_GBPD\an3to4_amortizedDA_GBCDtris.mat', 'data')
% load('H:\_190730\pseudo_GBPD\an3to4_amortizedDA_maskTouchFS.mat')
% data_all = [data_all; data];
% mask_trouch_fs_tri_all = [mask_trouch_fs_tri_all; mask_touch_fs_tri];

load('H:\_190730\pseudo_GBPD\An4new6_fixedOrigin_smooth_amortizedDA_GBCDtris.mat', 'data')
load('H:\_190730\pseudo_GBPD\An4new6_fixedOrigin_smooth_amortizedDA_maskTouchFS.mat')
data_all = [data_all; data];
mask_trouch_fs_tri_all = [mask_trouch_fs_tri_all; mask_touch_fs_tri];

fname = 'H:\_190730\pseudo_GBPD\allExcept3-4_amortizedDA';

% mask = data_all(:, end) > 0;
% csvwrite([fname, '_inner_pos.csv'],data_all(mask, :));
% csvwrite([fname, '_inner_neg.csv'],data_all(~mask, :));
% csvwrite([fname, '_inner.csv'],data_all);

data_all_fs = data_all(mask_trouch_fs_tri_all, :);
mask = data_all_fs(:, end) > 0;
csvwrite([fname, '_touchSF_pos.csv'],data_all_fs(mask, :));
csvwrite([fname, '_touchSF_neg.csv'],data_all_fs(~mask, :));
csvwrite([fname, '_touchSF.csv'],data_all_fs);

data_all_no_fs = data_all(~mask_trouch_fs_tri_all, :);
mask = data_all_no_fs(:, end) > 0;
csvwrite([fname, '_noSF_pos.csv'],data_all_no_fs(mask, :));
csvwrite([fname, '_noSF_neg.csv'],data_all_no_fs(~mask, :));
csvwrite([fname, '_noSF.csv'],data_all_no_fs);

%% #################################### Write file names ####################################
fname = 'file_names.txt';

format = '%s\n';
dictory = '..\pseudo_GBPD\';
states = {'An4new6_fixedOrigin_smooth', 'all_states', 'an0', 'an1with2', 'an2to3', 'an3to4'};
properties = {'_amortizedDA_complete', '_amortizedDA_complete_pos', '_amortizedDA_complete_neg', ...
    '_amortizedDA_incomplete', '_amortizedDA_incomplete_pos', '_amortizedDA_incomplete_neg',...
    '_amortizedDA_inner','_amortizedDA_inner_pos','_amortizedDA_inner_neg'};

fileID = fopen(fname,'w');
for i = 1:length(states)
    for j = 1:length(properties)
        name = [dictory, states{i}, properties{j}, '.csv'];
        fprintf(fileID, format, name);
    end
end
fclose('all');

%% #################################### Modify data in D3D ####################################
w_fa_diff = zeros(num_tris, 1);
w_fa_absdiff = zeros(num_tris, 1);
w_fa = zeros(num_tris, 1);
w_fa_diff_ratio = zeros(num_tris, 1);
w_fa_absdiff_ratio = zeros(num_tris, 1);
w_fitg_curv_diff = zeros(num_tris, 1);
w_figt_curv_absdiff = zeros(num_tris, 1);
w_fa_diff(idxes) = tri_amortized_a_diff;
w_fa_absdiff(idxes) = fa_absdiff;
w_fa(idxes) = tri_fa;
w_fa_diff_ratio(idxes) =  fa_diff_ratio;
w_fa_absdiff_ratio(idxes) =  fa_absdiff_ratio;
w_fitg_curv_diff(idxes) = fitg_curv_diff;
w_figt_curv_absdiff(idxes) = figt_curv_absdiff;

% h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/DAResidentFace', w_fa_diff');
% h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/DMsfResidentFace', w_fitg_curv_diff');
% h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/AreaResidentFace', w_fa');
% h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/DAResidentFaceRatio_noAreaLess20', w_fa_diff_ratio');
% h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/absDAResidentFaceRatio_noAreaLess20', w_fa_absdiff_ratio');
% h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/absDAResidentFace', w_fa_absdiff');
% h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/absDMsfResidentFace', w_figt_curv_absdiff');



%% #################################### Checks ####################################
file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
load('/Volumes/XIAOTING/Ni/190425_Hsmooth_geo_topo_an5crop2.mat', ...
    'face_area_an4', 'face_area_diff', 'face_itg_abscurv_an4', 'face_itg_abscurv_diff', 'tracked_uniqueface_an4');
mask = (abs(face_area_an4 + face_area_amortized_diff) > 10e-3 & face_area_an4 > 0);
face_area_an4 = face_area_an4(mask);
face_area_amortized_diff = face_area_amortized_diff(mask);
face_itg_abscurv_an4 = face_itg_abscurv_an4(mask);
face_itg_abscurv_diff = face_itg_abscurv_diff(mask);
obj_faces_an4 = obj_faces_an4(mask, :);


% -------- make keys --------
fl_key = cell(size(obj_faces_an4, 1), 1);
for i = 1:size(obj_faces_an4)
    %fl_key{i} = mat2str(tracked_uniqueface_an4(i, :));
    fl_key{i} = sprintf('%05d%05d', obj_faces_an4(i, 1), obj_faces_an4(i, 2));
end
face_area_an4_dict = containers.Map(fl_key, face_area_an4);
face_area_diff_dict = containers.Map(fl_key, face_area_amortized_diff);
face_area_diff_ratio = face_area_amortized_diff ./ face_area_an4;
face_area_diff_ratio_dict = containers.Map(fl_key, face_area_diff_ratio);

fl_d3d = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
face_area_diff_d3d = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/DAResidentFace')';
face_area_an4_d3d = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/AreaResidentFace');
face_area_diff_ratio_d3d = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/DAResidentFaceRatio')';
abs_face_area_diff_d3d = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/absDAResidentFace')';
abs_face_area_diff_ratio_d3d = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/absDAResidentFaceRatio')';

mask = face_area_an4_d3d > 0;
fl_d3d = fl_d3d(mask, :);
fl_d3d = sort(fl_d3d, 2);
face_area_diff_d3d = face_area_diff_d3d(mask);
face_area_an4_d3d = face_area_an4_d3d(mask);
face_area_diff_ratio_d3d = face_area_diff_ratio_d3d(mask);
abs_face_area_diff_d3d = abs_face_area_diff_d3d(mask);
abs_face_area_diff_ratio_d3d = abs_face_area_diff_ratio_d3d(mask);

sum(abs_face_area_diff_d3d - abs(face_area_diff_d3d))
sum(abs_face_area_diff_ratio_d3d - abs(face_area_diff_ratio_d3d))

data = [face_area_diff_ratio, face_area_an4, face_area_amortized_diff];
data = sortrows(data);
data(:,2) - data(:,3);

%%

rng('shuffle')
idx = randi(size(fl_d3d, 1));
disp('------------------------------')
disp(['facelabel = [', num2str(fl_d3d(idx, 1)), ', ', num2str(fl_d3d(idx, 2)), ']'])
disp(['A_tri = ', num2str(face_area_an4_d3d(idx)), ';    DA_tri = ', num2str(face_area_diff_d3d(idx)), ...
    ';  DAratio_tri = ', num2str(face_area_diff_ratio_d3d(idx))])
idx_fl = sprintf('%05d%05d', fl_d3d(idx, 1), fl_d3d(idx, 2));
disp(['A_face = ', num2str(face_area_an4_dict(idx_fl)), ';    DA_face = ', num2str(face_area_diff_dict(idx_fl)), ...
    ';  DAratio_face = ', num2str(face_area_diff_ratio_dict(idx_fl))])


%%
eps = 10e-3;
for i = 1:10000
    idx = randi(size(fl_d3d, 1));
    idx_fl = sprintf('%05d%05d', fl_d3d(idx, 1), fl_d3d(idx, 2));
    
    bad_value_1 = (face_area_an4_d3d(idx) - face_area_an4_dict(idx_fl)) > eps;    
    bad_value_2 = (face_area_diff_d3d(idx) - face_area_diff_dict(idx_fl)) > eps;    
    bad_value_3 = (face_area_diff_ratio_d3d(idx) - face_area_diff_ratio_dict(idx_fl)) > eps;    
    bad_value = bad_value_1 || bad_value_2 || bad_value_3;
    
    if bad_value
        warning(idx_fl);
    end
end


















