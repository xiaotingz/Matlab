% % ################################# Prepare data #################################
% tmp = readtable('/Volumes/XIAOTING/Ni/an0-an4/All_no_multi.xlsx','Sheet','An0-An1');
% look_up_table = table2array(tmp(:,1:2));
% save('look_up_table_an0_an1','look_up_table')
% 
% tmp = readtable('/Volumes/XIAOTING/Ni/an0-an4/All_no_multi.xlsx','Sheet','An1-An2');
% look_up_table = table2array(tmp(:,1:2));
% save('look_up_table_an1_an2','look_up_table')
% 
% tmp = readtable('/Volumes/XIAOTING/Ni/an0-an4/All_no_multi.xlsx','Sheet','An2-An3');
% look_up_table = table2array(tmp(:,1:2));
% save('look_up_table_an2_an3','look_up_table')
% 
% tmp = readtable('/Volumes/XIAOTING/Ni/an0-an4/All_no_multi.xlsx','Sheet','An3-An4');
% look_up_table = table2array(tmp(:,1:2));
% save('look_up_table_an3_an4','look_up_table')

file_an0 = '/Volumes/XIAOTING/Ni/an0-an4/An0new6.dream3d';
file_an1 = '/Volumes/XIAOTING/Ni/an0-an4/An1new6.dream3d';
file_an2 = '/Volumes/XIAOTING/Ni/an0-an4/An2new6.dream3d';
file_an3 = '/Volumes/XIAOTING/Ni/an0-an4/An3new6.dream3d';
% file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixedOrigin_smooth.dream3d';
file_an4 = '/Volumes/XIAOTING/Ni/an0-an4/An4new6.dream3d';
% file_an5_crop = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth_forParaview.dream3d';
file_an5 = '/Volumes/XIAOTING/Ni/an0-an4/An5new6.dream3d';
load('look_up_table_an0_an1');
look_up_table_an0_an1 = look_up_table;
load('look_up_table_an1_an2');
look_up_table_an1_an2 = look_up_table;
load('look_up_table_an2_an3');
look_up_table_an2_an3 = look_up_table;
load('look_up_table_an3_an4');
look_up_table_an3_an4 = look_up_table;
load('look_up_table_an4_an5crop')
look_up_table_an4_an5crop = look_up_table;
load('look_up_table_an4_an5')
look_up_table_an4_an5 = look_up_table;


dim_an0 = h5read(file_an0, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS')';
dim_an1 = h5read(file_an1, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS')';
dim_an2 = h5read(file_an2, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS')';
dim_an3 = h5read(file_an3, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS')';
dim_an4 = h5read(file_an4, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS')';
dim_an5 = h5read(file_an5, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS')';
steps = h5read(file_an0, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/SPACING')';

%% ################################# Rewrite file orign #################################
% h5write(file_an0, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', [0;0;0]);
% h5write(file_an1, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', [0;0;0]);
% h5write(file_an2, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', [0;0;0]);
% h5write(file_an3, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', [0;0;0]);
% h5write(file_an4, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', [0;0;0]);
% h5write(file_an5, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', [0;0;0]);


%% ################################# Fix shifts #################################
% avg_centroid_diff_an5_an4 = calcOriginShiftFromTrackedGrains(file_an5, file_an4, ...
%     [look_up_table_an4_an5(:,2), look_up_table_an4_an5(:,1)]);
% disp(avg_centroid_diff_an5_an4 ./ steps)
% shift_an5_an4 = calcVolumeDisplacement(file_an5, file_an4, ...
%     [look_up_table_an4_an5(:,2), look_up_table_an4_an5(:,1)], 20);


% avg_centroid_diff_an1_an0 = calcOriginShiftFromTrackedGrains(file_an1, file_an0, ...
%     [look_up_table_an0_an1(:,2), look_up_table_an0_an1(:,1)]);
% disp(avg_centroid_diff_an1_an0 ./ steps)
% shift_an1_an0 = calcVolumeDisplacement(file_an1, file_an0, ...
%     [look_up_table_an0_an1(:,2), look_up_table_an0_an1(:,1)], 20);
[volume_an1to0, volume_an0] = adjustSearchRange(shift_an1_an0, dim_an1, dim_an0);
disp('volume_d3d_an0')
disp(volume_an0 - 1)
disp('volume_d3d_an1to0')
disp(volume_an1to0 - 1)

% avg_centroid_diff_an2_an1 = calcOriginShiftFromTrackedGrains(file_an2, file_an1, ...
%     [look_up_table_an1_an2(:,2), look_up_table_an1_an2(:,1)]);
% disp(avg_centroid_diff_an2_an1 ./ steps)
% shift_an2_an1 = calcVolumeDisplacement(file_an2, file_an1, ...
%     [look_up_table_an1_an2(:,2), look_up_table_an1_an2(:,1)], 20);
[volume_an2to1, volume_an1with2] = adjustSearchRange(shift_an2_an1, dim_an2, dim_an1);
disp('volume_an1with2')
disp(volume_an1with2 - 1)
disp('volume_d3d_an2to1')
disp(volume_an2to1 - 1)

% avg_centroid_diff_an2_an3 = calcOriginShiftFromTrackedGrains(file_an2, file_an3, look_up_table_an2_an3);
% disp(avg_centroid_diff_an2_an3 ./ steps)
% shift_an2_an3 = calcVolumeDisplacement(file_an2, file_an3, look_up_table_an2_an3, 20);
[volume_an2to3, volume_an3with2] = adjustSearchRange(shift_an2_an3, dim_an2, dim_an3);
disp('volume_d3d_an2to3')
disp(volume_an2to3 - 1)
disp('volume_an3with2')
disp(volume_an3with2 - 1)
% 
% avg_centroid_diff_an3_an4 = calcOriginShiftFromTrackedGrains(file_an3, file_an4, look_up_table_an3_an4);
% disp(avg_centroid_diff_an3_an4 ./ steps)
% shift_an3_an4 = calcVolumeDisplacement(file_an3, file_an4, look_up_table_an3_an4, 20);
[volume_an3to4, volume_an4with3] = adjustSearchRange(shift_an3_an4, dim_an3, dim_an4);
disp('volume_d3d_an3to4')
disp(volume_an3to4 - 1)
disp('volume_an4with3')
disp(volume_an4with3 - 1)


%% ################################# Prepare look_up_table_crop #################################
file_an0 = '/Volumes/XIAOTING2/Ni_an0-an4/an0.dream3d';
file_an1to0 = '/Volumes/XIAOTING2/Ni_an0-an4/an1to0.dream3d';
file_an1with2 = '/Volumes/XIAOTING2/Ni_an0-an4/an1with2.dream3d';
file_an2to1 = '/Volumes/XIAOTING2/Ni_an0-an4/an2to1.dream3d';
file_an2to3 = '/Volumes/XIAOTING2/Ni_an0-an4/an2to3.dream3d';
file_an3with2 = '/Volumes/XIAOTING2/Ni_an0-an4/an3with2.dream3d';
file_an3to4 = '/Volumes/XIAOTING2/Ni_an0-an4/an3to4.dream3d';
file_an4with3 = '/Volumes/XIAOTING2/Ni_an0-an4/an4with3.dream3d';

% look_up_table_an0_an1crop = updateLookUpTableAfterCrop(file_an1, file_an1to0, 2, volume_an1to0, look_up_table_an0_an1);
% look_up_table = look_up_table_an0_an1crop;
% save('look_up_table_an0_an1crop.mat', 'look_up_table')

% look_up_table_an1crop_an2 = updateLookUpTableAfterCrop(file_an1, file_an1with2, 1, volume_an1with2, look_up_table_an1_an2);
% look_up_table_an1crop_an2crop = updateLookUpTableAfterCrop(file_an2, file_an2to1, 2, volume_an2to1, look_up_table_an1crop_an2);
% look_up_table = look_up_table_an1crop_an2crop;
% save('look_up_table_an1crop_an2crop.mat', 'look_up_table')


% look_up_table_an2crop_an3 = updateLookUpTableAfterCrop(file_an2, file_an2to3, 1, volume_an2to3, look_up_table_an2_an3);
% look_up_table_an2crop_an3crop = updateLookUpTableAfterCrop(file_an3, file_an3with2, 2, volume_an3with2, look_up_table_an2crop_an3);
% look_up_table = look_up_table_an2crop_an3crop;
% save('look_up_table_an2crop_an3crop.mat', 'look_up_table')

% look_up_table_an3crop_an4 = updateLookUpTableAfterCrop(file_an3, file_an3to4, 1, volume_an3to4, look_up_table_an3_an4);
% look_up_table_an3crop_an4crop = updateLookUpTableAfterCrop(file_an4, file_an4with3, 2, volume_an4with3, look_up_table_an3crop_an4);
% look_up_table = look_up_table_an3crop_an4crop;
% save('look_up_table_an3crop_an4crop.mat', 'look_up_table')

% save('look_up_tables_crop', 'look_up_table_an0_an1crop', 'look_up_table_an1crop_an2crop', ...
%     'look_up_table_an2crop_an3crop', 'look_up_table_an3crop_an4crop')



% ################################# distance from misorientations #################################
dg_obj = zeros(3, 3, 6);
dg_obj(:,:,1) = AAToG(60, [1, 1, 1]);     % sigma3
dg_obj(:,:,2) = AAToG(36.86, [1, 0, 0]);  % sigma5
dg_obj(:,:,3) = AAToG(38.21, [1, 1, 1]);  % sigma7
dg_obj(:,:,4) = AAToG(38.94, [1, 1, 0]);  % sigma9
dg_obj(:,:,5) = AAToG(31.59, [1, 1, 0]);  % sigma27a
dg_obj(:,:,6) = AAToG(35.43, [2, 1, 0]);  % sigma27b

% file = file_an0;
% load('/Volumes/XIAOTING2/Ni_an0-an4/geo_topo_an0_an1crop_full.mat', 'tracked_uniqueface_an4')
% dists_an0 = calcDistFromMisorientation(file, tracked_uniqueface_an4, dg_obj);
% 
% file = file_an1with2;
% load('/Volumes/XIAOTING2/Ni_an0-an4/geo_topo_an1crop_an2crop_full.mat', 'tracked_uniqueface_an4')
% dists_an1 = calcDistFromMisorientation(file, tracked_uniqueface_an4, dg_obj);
% 
% file = file_an2to3;
% load('/Volumes/XIAOTING2/Ni_an0-an4/geo_topo_an2crop_an3crop_full.mat', 'tracked_uniqueface_an4')
% dists_an2 = calcDistFromMisorientation(file, tracked_uniqueface_an4, dg_obj);
% 
% file = file_an3to4;
% load('/Volumes/XIAOTING2/Ni_an0-an4/geo_topo_an3crop_an4crop_full.mat', 'tracked_uniqueface_an4')
% dists_an3 = calcDistFromMisorientation(file, tracked_uniqueface_an4, dg_obj);
% 
% save('dists_all_states_full.mat', 'dists_an0', 'dists_an1', 'dists_an2', 'dists_an3');

file = file_an0;
load('/Volumes/XIAOTING2/Ni_an0-an4/geo_topo_an0_an1crop_inner.mat', 'tracked_uniqueface_an4')
dists_an0 = calcDistFromMisorientation(file, tracked_uniqueface_an4, dg_obj);

file = file_an1with2;
load('/Volumes/XIAOTING2/Ni_an0-an4/geo_topo_an1crop_an2crop_inner.mat', 'tracked_uniqueface_an4')
dists_an1 = calcDistFromMisorientation(file, tracked_uniqueface_an4, dg_obj);

file = file_an2to3;
load('/Volumes/XIAOTING2/Ni_an0-an4/geo_topo_an2crop_an3crop_inner.mat', 'tracked_uniqueface_an4')
dists_an2 = calcDistFromMisorientation(file, tracked_uniqueface_an4, dg_obj);

file = file_an3to4;
load('/Volumes/XIAOTING2/Ni_an0-an4/geo_topo_an3crop_an4crop_inner.mat', 'tracked_uniqueface_an4')
dists_an3 = calcDistFromMisorientation(file, tracked_uniqueface_an4, dg_obj);

save('dists_all_states_inner.mat', 'dists_an0', 'dists_an1', 'dists_an2', 'dists_an3');




%% ################################# combine geo_topo data #################################

face_area_init = [face_area_an4_an0; face_area_an4_an1; face_area_an4_an2; ...
    face_area_an4_an3; face_area_an4_an4];
face_area_diff = [face_area_diff_an0; face_area_diff_an1; face_area_diff_an2;...
    face_area_diff_an3; face_area_diff_an4];
face_itg_abscurv_init= [face_itg_abscurv_an4_an0; face_itg_abscurv_an4_an1; ...
    face_itg_abscurv_an4_an2; face_itg_abscurv_an4_an3; face_itg_abscurv_an4_an4];
face_itg_abscurv_diff = [face_itg_abscurv_diff_an0; face_itg_abscurv_diff_an1; ...
    face_itg_abscurv_diff_an2; face_itg_abscurv_diff_an3; face_itg_abscurv_diff_an4];


















