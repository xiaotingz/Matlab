% ##########################################################################
% * Notes
%     - Misorientaion distance is calculated from misorientation angle. 
%     - Distance to coherent twin is calculated as defined by Moraweic.
%         1. Identify (inner) complete faces 
%         2. Identify broken faces
%         3.1. Identify twins, from D3D
%         3.2. Myself distance from twins
%         4. Identify distance from directions
%         5. Write txt file
% ##########################################################################

% load('/Users/xiaotingzhong/Documents/Matlab/Grain Tracking/data/181108_rfvec.mat', 'rfvecs_an4', 'rfvecs_an5');
% file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% file_an5 = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth.dream3d';
file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin3_Hsmooth.dream3d';
file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4_Hsmooth.dream3d';
load('look_up_table_an4_an5crop.mat')

% ------------- get faces to track -------------
% [tracked_uniqueface_an4, tracked_uniqueface_an5] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_complete_faces');
% -------------
% load('/Volumes/XIAOTING/Ni/working/190621_tracked_faces_full.mat')
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190624_faces_full.mat')
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190425_Hsmooth_geo_topo_an5crop2.mat', ...
        'tracked_uniqueface_an4')
inner = ismember(sort(tracked_uniqueface_an4_full, 2), sort(tracked_uniqueface_an4, 2), 'rows');
% -------------

eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;


% % #################################### 1. Identify (inner) complete faces ####################################








% % #################################### 2. Identify broken faces ####################################
% is_onepiece = ones(size(tracked_uniqueface_an4, 1), 1);
% is_onepiece(face_piecewise) = 0;


% #################################### 3.1. Identify twins, from D3D ####################################
[is_twin_an4, incoherence_an4] = calcTwinsFromD3D(file_an4, faces_an4, eps_curv, eps_area, eps_min_ang);
[is_twin_an5, incoherence_an5] = calcTwinsFromD3D(file_an5, faces_an5, eps_curv, eps_area, eps_min_ang);

%%
% % #################################### 3.2. Myself distance from twins ####################################
% ----- By Rodrigues vector -----
% rfvec_twin = [1,1,1]/norm([1,1,1]) * tand(60/2);
% % """ Note 'AvgEAs' | 'AvgEulerAngles' """
% [rfvecs_an4]  = getFaceRFvecs(file_an4, tracked_uniqueface_an4);
% [rfvecs_an5]  = getFaceRFvecs(file_an5, tracked_uniqueface_an5);
% dist_twin_an4 = vecnorm(rfvecs_an4 - rfvec_twin, 2, 2);
% dist_twin_an5 = vecnorm(rfvecs_an5 - rfvec_twin, 2, 2);

% ----- By misorientation matrix -----
dg_obj = zeros(3, 3, 6);
dg_obj(:,:,1) = AAToG(60, [1, 1, 1]);     % sigma3
dg_obj(:,:,2) = AAToG(36.86, [1, 0, 0]);  % sigma5
dg_obj(:,:,3) = AAToG(38.21, [1, 1, 1]);  % sigma7
dg_obj(:,:,4) = AAToG(38.94, [1, 1, 0]);  % sigma9
dg_obj(:,:,5) = AAToG(31.59, [1, 1, 0]);  % sigma27a
dg_obj(:,:,6) = AAToG(35.43, [2, 1, 0]);  % sigma27b

dists_miso = calcDistFromMisorientation(file_an4, faces_an4, dg_obj);


%% #################################### 4. Identify distance from directions ####################################
target_normals = [1,1,1;1,1,0;1,0,0;1,1,2];
target_normals = target_normals ./ vecnorm(target_normals, 2, 2);
[dists_norm_an4, std_norm_an4] = calcDistanceFromPlaneNormals(file_an4, faces_an4, target_normals, eps_curv, eps_area, eps_min_ang);


%% #################################### 5. Write txt file ####################################
% not_twin_an4 = ~istwin_an4;
% not_twin_an5 = ~istwin_an5;
% dist_twin_an4 = dists_an4_full(:,1);

fileID = fopen('190624_features_otherinfo.txt','w');
% fprintf(fileID,'%s, %s, %s, %s\n','is_onepiece', 'not_twin', 'dist_twin', 'weighteddist_ctwin');
fprintf(fileID,'%s, %s, %s, %s\n','not_twin_an4', 'dist_twin_an4', 'incoherence_an4', 'mask_inner');
for i = 1:length(not_twin_an4)
    if mask_good_face(i)
        fprintf(fileID, '%6d, %6.3f, %6.3f, %6d\n', ...
            not_twin_an4(i), dist_twin_an4(i), incoherence_an4(i), inner(i)) ;
    end
end
fclose(fileID);

