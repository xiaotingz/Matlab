%% #################################### Data Prepare ####################################
file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
load('/Volumes/XIAOTING/Ni/190425_Hsmooth_geo_topo_an5crop2.mat', ...
    'face_area_an4', 'face_area_diff', 'face_itg_abscurv_an4', 'face_itg_abscurv_diff', 'tracked_uniqueface_an4');

%  ------------------------------------- read and clean data -------------------------------------
fl = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
curv = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
min_ang = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
area = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
eps_area = 7;
eps_curv = 1;
eps_min_ang = 10;
num_tris = length(curv);

fl = sort(fl, 2);
idxes = (1:length(fl))';
mask = all(fl >= 0, 2);
fl = fl(mask, :);
curv = curv(mask);
min_ang = min_ang(mask);
area = area(mask);
idxes = idxes(mask);

valid_grainids = unique(tracked_uniqueface_an4);
mask = all(fl>0, 2) & abs(curv) < eps_curv & area < eps_area & min_ang > eps_min_ang ...
        & all(ismember(fl, valid_grainids), 2);
% data = data(mask, :);
fl = fl(mask, :);
area = area(mask, :);
curv = curv(mask);
idxes = idxes(mask);

%% #################################### Assign Resident Face Values to Triangles ####################################

% %  ------------------------------------- create dictionaries -------------------------------------
% fl_keys = cell(length(tracked_uniqueface_an4),1);
% for i = 1:length(tracked_uniqueface_an4)
%     fl_keys{i} = mat2str(tracked_uniqueface_an4(i, :));
% end
% fa_dict = containers.Map(fl_keys,face_area_an4);
% fa_diff_dict = containers.Map(fl_keys,face_area_diff);
% fitg_curv_dict = containers.Map(fl_keys,face_itg_abscurv_an4);
% fitg_curv_diff_dict = containers.Map(fl_keys,face_itg_abscurv_diff);


%  ------------------------------------- assign face values to individual triangles -------------------------------------
fa = zeros(size(fl, 1), 1) * NaN;
fa_diff = zeros(size(fl, 1), 1) * NaN;
fitg_curv = zeros(size(fl, 1), 1) * NaN;
fitg_curv_diff = zeros(size(fl, 1), 1) * NaN;

% for i = 1:length(fl)
%     fl_key = mat2str(fl(i, :));
%     if isKey(fa_dict, fl_key)
%         fa(i) = fa_dict(fl_key);
%         fa_diff(i) = fa_diff_dict(fl_key);
%         fitg_curv(i) = fitg_curv_dict(fl_key);
%         fitg_curv_diff(i) = fitg_curv_diff_dict(fl_key);
%     end
% end
for i = 1:length(tracked_uniqueface_an4)
    mask = (fl(:, 1) == tracked_uniqueface_an4(i, 1) & fl(:, 2) == tracked_uniqueface_an4(i, 2));
    fa(mask) = face_area_an4(i);
    fa_diff(mask) = face_area_diff(i);
    fitg_curv(mask) = face_itg_abscurv_an4(i);
    fitg_curv_diff(mask) = face_itg_abscurv_diff(i);
end


mask = ~isnan(fa);
fa = fa(mask);
fa_diff = fa_diff(mask);
fa_absdiff = abs(fa_diff);
fa_diff_ratio = fa_diff ./ fa;
fa_absdiff_ratio = fa_absdiff ./ fa;
fitg_curv = fitg_curv(mask);
fitg_curv_diff = fitg_curv_diff(mask);
figt_curv_absdiff = abs(fitg_curv_diff);
idxes = idxes(mask);

% data = data(mask, :);
curv = curv(mask);


%% #################################### Write Seudo GBCDtris file ####################################
% fa_diff_id = fopen('an4_fa_diff.txt', 'w');
% format = '10%6.3f';
% fprintf(fa_diff_id,format, data())
% data(:, end) = fa_diff;
% csvwrite('an4_fa_diff.csv',data);

data(:, end) = fa_absdiff;
csvwrite('an4_fa_absdiff.csv',data);



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

h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/DAResidentFace', w_fa_diff');
h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/DMsfResidentFace', w_fitg_curv_diff');
h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/AreaResidentFace', w_fa');
h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/DAResidentFaceRatio', w_fa_diff_ratio');
h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/absDAResidentFaceRatio', w_fa_absdiff_ratio');
h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/absDAResidentFace', w_fa_absdiff');
h5write(file_an4, '/DataContainers/TriangleDataContainer/FaceData/absDMsfResidentFace', w_figt_curv_absdiff');






