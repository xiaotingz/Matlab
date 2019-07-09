function look_up_table = updateLookUpTableAfterCrop(file_full, file_crop, crop_idx, volume, look_up_table)
% ##########################################################################
% * Input
%     - crop_idx = 1 or 2,
%             Indicating which file was cropped, and which column of 
%             look_up_table to change.
%     - volume = [3, 2] 
%             The range of cropped volume. 
%             Returned by adjustSearchRange.m. Or see prepareAllStates.m for example.
% * Output
%     - look_up_table = [n, 2]  
%             The updated look_up_table
% * Objective of the function
%     - The different anneal states need to be cropped to obtain pairwise common volumes.
%       Cropping renumber grains. This function updates the look_up_table so we
%       do not lose track of the grains.
% ##########################################################################
% ----------------------- load debug data -----------------------
% file_full = '/Volumes/XIAOTING/Ni/An5new6_Hsmooth.dream3d';
% file_crop = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth.dream3d';
% file_full = '/Volumes/XIAOTING/Ni/an0-an4/An1new6.dream3d';
% file_crop = '/Volumes/XIAOTING2/Ni_an0-an4/an1to0.dream3d';
% volume = volume_an1to0;
% crop_idx = 2;
% ---------------------------------------------------------------
%% ##### Construct a Map (Dictionary) that converts the ids #####
% """ 
% check if voxels the cropped volume has a 1-to-1 corresp with the 
% intended volume in an5. If yes, D3D crop process correct.
% """
voxel_gid_full_raw = squeeze(h5read(file_full, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));
voxel_gid_crop = squeeze(h5read(file_crop, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));
% voxel_gid_full = voxel_gid_full_raw(4:398, :, 4:87);
voxel_gid_full = voxel_gid_full_raw(volume(1,1):volume(1,2), ...
    volume(2,1):volume(2,2), volume(3,1):volume(3,2));

map_old_to_crop = containers.Map('KeyType','int32','ValueType','int32');
seen = [];
for i = 1:size(voxel_gid_crop, 1)
    for j = 1:size(voxel_gid_crop, 2)
        for k = 1:size(voxel_gid_crop, 3)
            id_old = voxel_gid_full(i, j, k);
            id_new = voxel_gid_crop(i, j, k);
            
            if ~ismember(id_old, seen)
                map_old_to_crop(id_old) = id_new;
                seen = [seen; id_old];
            else
                if map_old_to_crop(id_old) ~= id_new
                    warning('wrong!');
                    break
                end
            end
        end
    end
end

% ############### Construct a new look_up_table between an4 & an5_crop ###############
% ----- first remove the tracked grains that got cropped away -----
% """ voxel_gid_full contain id=0 but look_up_table don't. So it's fine. """
mask = ismember(look_up_table(:, crop_idx), unique(voxel_gid_full));
look_up_table = look_up_table(mask, :);

% ----- then convert an5_id, from old_full to new_crop -----
for i = 1:length(look_up_table)
    look_up_table(i, crop_idx) = map_old_to_crop(look_up_table(i, crop_idx));
end

% ----- sort look_up_table -----
look_up_table = sortrows(look_up_table, 2);

end

% % ------- check function -------
% file_an5 = '/Volumes/XIAOTING/Ni/an0-an4/An5new6.dream3d';
% file_an5crop = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth.dream3d';
% load('look_up_table_an4_an5.mat')
% [volume_an5to4, ~] = adjustSearchRange([3,0,3], dim_an5, dim_an4);
% look_up_table_an4_an5crop = updateLookUpTableAfterCrop(file_an5, file_an5crop, 2, volume_an5to4, look_up_table);
% load('look_up_table_an4_an5crop.mat')
% sum(look_up_table == look_up_table_an4_an5crop)
% % ------------------------------


% %% ############### Double Check the Map ###############
% voxel_gid_cropfromfull = -1 * ones(size(voxel_gid_crop));
% for i = 1:size(voxel_gid_crop, 1)
%     for j = 1:size(voxel_gid_crop, 2)
%         for k = 1:size(voxel_gid_crop, 3)
%             voxel_gid_cropfromfull(i, j, k) = map_old_to_crop(voxel_gid_full(i, j, k));
%         end
%     end
% end
% 
% tmp = sum(voxel_gid_cropfromfull ~= voxel_gid_crop);
% sum(tmp(:))
% 
% %% ############### Check in paraview ###############
% load('/Volumes/XIAOTING/Ni/working/190425_Hsmooth_geo_topo_an5crop2.mat', ...
%     'tracked_uniqueface_an4', 'tracked_uniqueface_an5');
% file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% file_an5 = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth.dream3d';
% load('look_up_table_an4_an5crop.mat')
% 
% rng('shuffle')
% idx = randi(length(look_up_table));
% 
% 
% % disp(look_up_table(idx, :));
% disp('-------------- ')
% dispFacePairInfo(file_an4, file_an5, tracked_uniqueface_an4, tracked_uniqueface_an5, idx)





