%% ##### Construct a Map (Dictionary) that converts the ids #####
file_full = '/Volumes/XIAOTING/Ni/An5new6_Hsmooth.dream3d';
file_crop = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth.dream3d';

voxel_gid_full_raw = squeeze(h5read(file_full, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));
voxel_gid_crop = squeeze(h5read(file_crop, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));
voxel_gid_full = voxel_gid_full_raw(4:398, :, 4:87);

map = containers.Map('KeyType','int32','ValueType','int32');
seen = [];
for i = 1:size(voxel_gid_crop, 1)
    for j = 1:size(voxel_gid_crop, 2)
        for k = 1:size(voxel_gid_crop, 3)
            id_old = voxel_gid_full(i, j, k);
            id_new = voxel_gid_crop(i, j, k);
            
            if ~ismember(id_old, seen)
                map(id_old) = id_new;
                seen = [seen; id_old];
            else
                if map(id_old) ~= id_new
                    warning('wrong!');
                    break
                end
            end
        end
    end
end

id_l_oldfull_r_newcrop =  [cell2mat(keys(map)); cell2mat(values(map))]';                
map_keyfull_valcrop = map;

%% ##### Double Check the Map #####
voxel_gid_cropfromfull = -1 * ones(size(voxel_gid_crop));
for i = 1:size(voxel_gid_crop, 1)
    for j = 1:size(voxel_gid_crop, 2)
        for k = 1:size(voxel_gid_crop, 3)
            voxel_gid_cropfromfull(i, j, k) = map(voxel_gid_full(i, j, k));
        end
    end
end

tmp = sum(voxel_gid_cropfromfull ~= voxel_gid_crop);
sum(tmp(:))

%% ##### Construct a new look_up_table between an4 & an5_crop #####
% ----- first remove the tracked an5 grains that got cropped away -----
mask = ismember(look_up_table(:, 2), id_1c_oldfull_2c_newcrop(:, 1));
look_up_table = look_up_table(mask, :);

% ----- then convert an5_id, from old_full to new_crop -----
for i = 1:length(look_up_table)
    look_up_table(i, 2) = map(look_up_table(i, 2));
end

% ----- sort look_up_table by id_an5 -----
look_up_table = sortrows(look_up_table, 2);


%% ##### Check in paraview #####
rng('shuffle')
idx = randi(length(look_up_table));


disp(look_up_table(idx, :))














