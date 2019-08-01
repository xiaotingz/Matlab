%% #################################### Read Data ####################################
% % ----------- inner faces -----------
% file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
% load('/Volumes/XIAOTING/Ni/working/190425_Hsmooth_geo_topo_an5crop2.mat', ...
%     'face_area_an4', 'face_area_diff', 'face_itg_abscurv_an4', 'face_itg_abscurv_diff', 'tracked_uniqueface_an4');
fname = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/GBCD_GBPD/pseudo_gbpd/Lsmooth_amortized';
% % ----------- full tracked faces -----------
% load('/Volumes/XIAOTING/Ni/working/190621_tracked_faces_full.mat')
% tracked_uniqueface_an4 = tracked_uniqueface_an4_full;
% tracked_uniqueface_an5 = tracked_uniqueface_an5_full;
% % ----------- all_states -----------
% data = an0GBCDtris;
% load('/Volumes/XIAOTING2/Ni_an0-an4/geo_topo_an0_an1crop_full.mat', 'face_area_an4', 'face_area_diff', 'tracked_uniqueface_an4');
% file_an4 = '/Volumes/XIAOTING2/Ni_an0-an4/an0.dream3d';
% fname = '/Volumes/XIAOTING2/Ni_an0-an4/fullTrack_an0_an1crop';
% % ----------- iron -----------
% file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Iron/iron_an0.dream3d';
% fname = '/Users/xiaotingzhong/Desktop/Datas/Iron/iron_amortized';

% """
% main.m to calculate face_itg_curv_an4
% """
face_area_an4 = face_itg_curv_an4(:,1);

face_area_diff_ratio = face_area_diff ./ face_area_an4;
mask = ((face_area_an4 + face_area_diff) > 20 & face_area_an4 > 20 & ... 
    face_area_diff_ratio > -0.9 & face_area_diff_ratio < 10);
face_area_an4 = face_area_an4(mask);
face_area_diff = face_area_diff(mask);
tracked_uniqueface_an4 = tracked_uniqueface_an4(mask, :);
face_area_diff_ratio = face_area_diff_ratio(mask);

% ---------------------------------- Amortize area change to each triangle ----------------------------------
face_num_tris_an4 = calcNumTrianglesOnFace(file_an4, tracked_uniqueface_an4, eps_curv, eps_area, eps_min_ang);
face_area_diff = face_area_diff./ face_num_tris_an4;
% -----------------------------------------------------------------------------------------------------------
%%

% %% #################################### Data Prepare, Hsmoosth ####################################
% % """
% % bad triangles in Hsmooth are found from eps_area, eps_curv and eps_min_ang
% % """
% fl = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
% curv = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
% min_ang = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
% area = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
% eps_area = 7;
% eps_curv = 1;
% eps_min_ang = 10;
% num_tris = length(curv);
% 
% fl = sort(fl, 2);
% idxes = (1:length(fl))';
% mask = all(fl >= 0, 2);
% fl = fl(mask, :);
% curv = curv(mask);
% min_ang = min_ang(mask);
% area = area(mask);
% 
% idxes = idxes(mask);
% valid_grainids = unique(tracked_uniqueface_an4);
% mask = all(fl>0, 2) & abs(curv) < eps_curv & area < eps_area & min_ang > eps_min_ang ...
%         & all(ismember(fl, valid_grainids), 2);
% data = data(mask, :);
% fl = fl(mask, :);
% area = area(mask, :);
% curv = curv(mask);
% idxes = idxes(mask);


%%% #################################### Data Prepare, Laplacian Smoosth ####################################
% """
% exclude triple line triangles
% """
tri_nodes = 1 + h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')';
node_types = h5read(file_an4,'/DataContainers/TriangleDataContainer/VertexData/NodeType')';
fl = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
fl = sort(fl, 2);
mask = all(fl >= 0, 2);
fl = fl(mask, :);
tri_nodes = tri_nodes(mask, :);

tri_node_types = node_types(tri_nodes);
mask_good_tris = (all(tri_node_types == 2, 2) & all(fl > 0, 2));

fl = fl(mask_good_tris, :);
data = data(mask_good_tris, :);


%% #################################### Assign Resident Face Values to Triangles ####################################
%  ------------------------------------- assign face values to individual triangles -------------------------------------
fa = zeros(size(fl, 1), 1) * NaN;
fa_diff = zeros(size(fl, 1), 1) * NaN;
fitg_curv = zeros(size(fl, 1), 1) * NaN;
fitg_curv_diff = zeros(size(fl, 1), 1) * NaN;

for i = 1:length(tracked_uniqueface_an4)
    mask = (fl(:, 1) == tracked_uniqueface_an4(i, 1) & fl(:, 2) == tracked_uniqueface_an4(i, 2));
    fa(mask) = face_area_an4(i);
    fa_diff(mask) = face_area_diff(i);
end


mask_trackedface = ~isnan(fa);
fa = fa(mask_trackedface);
fa_diff = fa_diff(mask_trackedface);
fa_absdiff = abs(fa_diff);
fa_diff_ratio = fa_diff ./ fa;
fa_absdiff_ratio = fa_absdiff ./ fa;
% idxes = idxes(mask_trackedface);

data = data(mask_trackedface, :);
% curv = curv(mask_trackedface);

save([fname, '_GBCDtris.mat'], 'data', 'fa', 'fa_diff', 'fa_absdiff', ...
    'fa_diff_ratio', 'fa_absdiff_ratio', 'face_area_an4', 'face_area_diff', ...
    'tracked_uniqueface_an4');


% #################################### Write Seudo GBCDtris file ####################################
% data = [data_an0; data_an1; data_an2; data_an3; data_an4];
% fa_diff = [fa_diff_an0; fa_diff_an1; fa_diff_an2; fa_diff_an3; fa_diff_an4];
% fa_diff_ratio = [fa_diff_ratio_an0; fa_diff_ratio_an1; fa_diff_ratio_an2; ...
%     fa_diff_ratio_an3; fa_diff_ratio_an4];
% fa_absdiff = [fa_absdiff_an0; fa_absdiff_an1; fa_absdiff_an2; fa_absdiff_an3; fa_absdiff_an4];
% fa_absdiff_ratio = [fa_absdiff_ratio_an0; fa_absdiff_ratio_an1; ...
%     fa_absdiff_ratio_an2; fa_absdiff_ratio_an3; fa_absdiff_ratio_an4];
% fname = 'fullTrack_all_states';

%%
% data = [data, fa_diff];
data(:, end) = fa_diff;
mask = fa_diff > 0;
csvwrite([fname,'_DA_pos.csv'],data(mask, :));
csvwrite([fname,'_DA_neg.csv'],data(~mask, :));
csvwrite([fname,'_DA.csv'],data);
% data(:, end) = fa_diff_ratio;
% csvwrite([fname,'_DAratio_pos.csv'],data(mask, :));
% csvwrite([fname,'_DAratio_neg.csv'],data(~mask, :));
% csvwrite([fname,'_DAratio.csv'],data);
% data(:, end) = fa_absdiff/100;
% csvwrite([fname,'_001absDA.csv'],data);
% data(:, end) = fa_absdiff_ratio;
% csvwrite([fname,'_absDAratio.csv'],data);


% data = [data, fa_diff/100];
% mask = fa_diff > 0;
% csvwrite('noExtreme_Lsmooth_001DA_pos.csv',data(mask, :));
% csvwrite('noExtreme_Lsmooth_001DA_neg.csv',data(~mask, :));
% csvwrite('noExtreme_Lsmooth_001DA.csv',data);
% data(:, end) = fa_diff_ratio;
% csvwrite('noExtreme_Lsmooth_DAratio_pos.csv',data(mask, :));
% csvwrite('noSExtreme_Lsmooth_DAratio_neg.csv',data(~mask, :));
% csvwrite('noExtreme_Lsmooth_DAratio.csv',data);
% data(:, end) = fa_absdiff/100;
% csvwrite('noExtreme_Lsmooth_001absDA.csv',data);
% data(:, end) = fa_absdiff_ratio;
% csvwrite('noExtreme_Lsmooth_absDAratio.csv',data);

% file = 'GBCDtris_noSmallFace_A.txt';
% fileID = fopen(file,'w');
% format = '%7.4f    %7.4f    %7.4f    %7.4f    %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f   %7.4f\n';
% for i = 1:length(data)
%     fprintf(fileID,format,data(i, 1),data(i, 2),data(i, 3),data(i, 4),data(i, 5),data(i, 6), ...
%         data(i, 7),data(i, 8),data(i, 9), data(i, 10), fa(i));
% end
% fclose('all');


% file = 'GBCDtris_001DA.txt';
% fileID = fopen(file,'w');
% format = '%7.4f    %7.4f    %7.4f    %7.4f    %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f   %7.4f\n';
% for i = 1:length(data)
%     fprintf(fileID,format,data(i, 1),data(i, 2),data(i, 3),data(i, 4),data(i, 5),data(i, 6), ...
%         data(i, 7),data(i, 8),data(i, 9),fa(i), fa_diff(i));
% end
% fclose('all');

% file = 'GBCDtris_absDA.txt';
% fileID = fopen(file,'w');
% format = '%7.4f    %7.4f    %7.4f    %7.4f    %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f   %7.4f\n';
% for i = 1:length(data)
%     fprintf(fileID,format,data(i, 1),data(i, 2),data(i, 3),data(i, 4),data(i, 5),data(i, 6), ...
%         data(i, 7),data(i, 8),data(i, 9),fa(i), abs(fa_diff(i)));
% end
% fclose('all');

% file = 'GBCDtris_DAratio.txt';
% fileID = fopen(file,'w');
% format = '%7.4f    %7.4f    %7.4f    %7.4f    %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f   %7.4f\n';
% for i = 1:length(data)
%     fprintf(fileID,format,data(i, 1),data(i, 2),data(i, 3),data(i, 4),data(i, 5),data(i, 6), ...
%         data(i, 7),data(i, 8),data(i, 9), data(i, 10), fa_diff_ratio(i));
% end
% fclose('all');
% 
% 
% file = 'GBCDtris_absDAratio.txt';
% fileID = fopen(file,'w');
% format = '%7.4f    %7.4f    %7.4f    %7.4f    %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f   %7.4f\n';
% for i = 1:length(data)
%     fprintf(fileID,format,data(i, 1),data(i, 2),data(i, 3),data(i, 4),data(i, 5),data(i, 6), ...
%         data(i, 7),data(i, 8),data(i, 9),data(i, 10), fa_absdiff_ratio(i));
% end
% fclose('all');


%% #################################### Modify data in D3D ####################################
w_fa_diff = zeros(num_tris, 1);
w_fa_absdiff = zeros(num_tris, 1);
w_fa = zeros(num_tris, 1);
w_fa_diff_ratio = zeros(num_tris, 1);
w_fa_absdiff_ratio = zeros(num_tris, 1);
w_fitg_curv_diff = zeros(num_tris, 1);
w_figt_curv_absdiff = zeros(num_tris, 1);
w_fa_diff(idxes) = fa_diff;
w_fa_absdiff(idxes) = fa_absdiff;
w_fa(idxes) = fa;
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
mask = (abs(face_area_an4 + face_area_diff) > 10e-3 & face_area_an4 > 0);
face_area_an4 = face_area_an4(mask);
face_area_diff = face_area_diff(mask);
face_itg_abscurv_an4 = face_itg_abscurv_an4(mask);
face_itg_abscurv_diff = face_itg_abscurv_diff(mask);
tracked_uniqueface_an4 = tracked_uniqueface_an4(mask, :);


% -------- make keys --------
fl_key = cell(size(tracked_uniqueface_an4, 1), 1);
for i = 1:size(tracked_uniqueface_an4)
    %fl_key{i} = mat2str(tracked_uniqueface_an4(i, :));
    fl_key{i} = sprintf('%05d%05d', tracked_uniqueface_an4(i, 1), tracked_uniqueface_an4(i, 2));
end
face_area_an4_dict = containers.Map(fl_key, face_area_an4);
face_area_diff_dict = containers.Map(fl_key, face_area_diff);
face_area_diff_ratio = face_area_diff ./ face_area_an4;
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

data = [face_area_diff_ratio, face_area_an4, face_area_diff];
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


















