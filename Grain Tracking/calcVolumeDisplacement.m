function shift = calcVolumeDisplacement(file_an4, file_an5, look_up_table)
% ##################################################################
% * Notes
%   - This function finds the best shift between two volumes by counting the
% #(overlapping voxels).
%   - Dependencies: adjustSearchRange.m
% ##################################################################
% ------------------------------------------------------
% file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin_centroid_smooth.dream3d');
% file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
% load('look_up_table_an4_an5.mat')
% ------------------------------------------------------
gid_an4 = squeeze(h5read(file_an4,'/DataContainers/ImageDataContainer/CellData/FeatureIds'));
gid_an5 = squeeze(h5read(file_an5,'/DataContainers/ImageDataContainer/CellData/FeatureIds'));
dim_an4 = size(gid_an4);
dim_an5 = size(gid_an5);

% ################################### Change FeatureID in an4 to Id of an5 ###################################
id_dic = containers.Map(look_up_table(:,1),look_up_table(:,2));
for i = 1:size(gid_an4, 1)
    for j = 1:size(gid_an4, 2)
        for k = 1:size(gid_an4, 3)
            if any(look_up_table(:,1) == gid_an4(i, j, k))
                gid_an4(i, j, k) = id_dic(gid_an4(i, j, k));
            else
                gid_an4(i, j, k) = -1;
            end
        end
    end
end

% """
% Find best position shift: coordinate descent. 
% """
shift = [0,0,0];
x_range_start = -2:5;
z_range_start = 0:20;
x_range_normal = -2:2;
z_range_normal = -3:3;
y_range = -2:2;
most_match = 0;
continue_opt = true;
first_search = true;

while continue_opt
    % ----- Get search range -----
    if first_search
        first_search = false;
        x_range = x_range_start;
        z_range = z_range_start;
    else
        x_range = x_range_normal;
        z_range = z_range_normal;
    end
    
    % ----- Record old shift value -----
    shift_old = shift;
    
    % ----- Search for a better shift: coordinate descent x, y, z -----
    for i = 1:length(x_range)
        shift_tmp = shift;
        shift_tmp(1) = x_range(i);
        [idx_an4, idx_an5] = adjustSearchRange(shift_tmp, dim_an4, dim_an5);
        mask = (gid_an4(idx_an4(1,1):idx_an4(1,2), idx_an4(2,1):idx_an4(2,2), idx_an4(3,1):idx_an4(3,2)) ...
                == gid_an5(idx_an5(1,1):idx_an5(1,2), idx_an5(2,1):idx_an5(2,2), idx_an5(3,1):idx_an5(3,2)));
        if sum(mask(:)) > most_match
            most_match = sum(mask(:));
            shift = shift_tmp;
        end
    end
    for j = 1:length(y_range)
        shift_tmp = shift;
        shift_tmp(2) = y_range(j);
        [idx_an4, idx_an5] = adjustSearchRange(shift_tmp, dim_an4, dim_an5);
        mask = (gid_an4(idx_an4(1,1):idx_an4(1,2), idx_an4(2,1):idx_an4(2,2), idx_an4(3,1):idx_an4(3,2)) ...
                == gid_an5(idx_an5(1,1):idx_an5(1,2), idx_an5(2,1):idx_an5(2,2), idx_an5(3,1):idx_an5(3,2)));
        if sum(mask(:)) > most_match
            most_match = sum(mask(:));
            shift = shift_tmp;
        end
    end
    for k = 1:length(z_range)
        shift_tmp = shift;
        shift_tmp(3) = z_range(k);
        [idx_an4, idx_an5] = adjustSearchRange(shift_tmp, dim_an4, dim_an5);
        mask = (gid_an4(idx_an4(1,1):idx_an4(1,2), idx_an4(2,1):idx_an4(2,2), idx_an4(3,1):idx_an4(3,2)) ...
                == gid_an5(idx_an5(1,1):idx_an5(1,2), idx_an5(2,1):idx_an5(2,2), idx_an5(3,1):idx_an5(3,2)));
        if sum(mask(:)) > most_match
            most_match = sum(mask(:));
            shift = shift_tmp;
        end
    end    
    
    % ----- Converge = the best shift didn't change -----
    if all(shift == shift_old)
        continue_opt = false;
    end
    
end

% disp('final')
% disp(shift)
end




