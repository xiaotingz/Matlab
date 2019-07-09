function [shift_best, ranges_start] = calcVolumeDisplacement(file_a, file_b, look_up_table, lim)
% ##################################################################
% * Notes
%   - This function finds the best shift between two volumes by counting the
% #(overlapping voxels).
%   - Input
%         - look_up_table = [id_file_a, id_file_b]
%         - file_a = smaller_volume, file_b = larger_volume is the best to
%         understand though not required. 
%                 Remember tp adjust look_up_table according to file_a and file_b.
%         - lim, a number from guessed from avg_tracked_inner_grains_centroid_diff
%                 For default, use 3. But it can be very large if field of
%                 view too different.
%   - Dependencies: adjustSearchRange.m
% ##################################################################
% ------------------------------------------------------
% file_a = ('/Users/xiaotingzhong/Desktop/Datas/Ni_a_5/An4new6_fixOrigin_centroid_smooth.dream3d');
% file_b = ('/Users/xiaotingzhong/Desktop/Datas/Ni_a_5/An5new6_smooth.dream3d');
% load('look_up_table_an4_an5.mat')
% file_a = file_an1;
% file_b = file_an2;
% x_range_start = ranges_start{1};
% y_range_start = ranges_start{2};
% z_range_start = ranges_start{3};
% look_up_table = look_up_table_an1_an2;
% ------------------------------------------------------
gid_a = squeeze(h5read(file_a,'/DataContainers/ImageDataContainer/CellData/FeatureIds'));
gid_b = squeeze(h5read(file_b,'/DataContainers/ImageDataContainer/CellData/FeatureIds'));
dim_a = size(gid_a);
dim_b = size(gid_b);

% ################################### Change FeatureID in an4 to Id of an5 ###################################
% """
% Find best position shift: coordinate descent. 
%         In the first search look all possible ranges. Later adjust locally. 
% """
id_dic = containers.Map(look_up_table(:,1),look_up_table(:,2));
for i = 1:size(gid_a, 1)
    for j = 1:size(gid_a, 2)
        for k = 1:size(gid_a, 3)
            if any(look_up_table(:,1) == gid_a(i, j, k))
                gid_a(i, j, k) = id_dic(gid_a(i, j, k));
            else
                gid_a(i, j, k) = -1;
            end
        end
    end
end


% x_range_start = -2:5;
% y_range_start = -3:3;
% z_range_start = 0:20;
% """
% Initial search range:
%       Note volume_a stable and volume_b moving, 
%       If dim_a > dim_b, should do more positive motion. 
%           volume_a    |--------------------------|
%           volume_b  - |----------------| +                        (mobile)
%       If dim_a < dim_b, should do more negative motion.
%           volume_a    |--------------------------|
%           volume_b  - |------------------------------------| +    (mobile)
% """
ranges_start = cell(3, 1);
for i = 1:3
    diff = abs(dim_b(i) - dim_a(i));
    if dim_a(i) > dim_b(i)
        ranges_start{i} = (-lim : diff + lim)';
    else
        ranges_start{i} = (-lim - diff :  lim)';
    end
end
x_range_normal = -3:3;
z_range_normal = -3:3;
y_range_normal = -3:3;

shift_best = [0,0,0];
most_match = 0;
continue_opt = true;
first_search = true;
while continue_opt
    % ----- Get search range -----
    if first_search
        first_search = false;
        x_range = ranges_start{1};
        y_range = ranges_start{2};
        z_range = ranges_start{3};
    else
        x_range = x_range_normal;
        y_range = y_range_normal;
        z_range = z_range_normal;
    end
    
    % ----- Record old shift value -----
    shift_old = shift_best;
    
    % ----- Search for a better shift: coordinate descent x, y, z -----
    for i = 1:length(x_range)
        shift_tmp = shift_best;
        shift_tmp(1) = x_range(i);
        [idx_a, idx_b] = adjustSearchRange(shift_tmp, dim_a, dim_b);
        mask = (gid_a(idx_a(1,1):idx_a(1,2), idx_a(2,1):idx_a(2,2), idx_a(3,1):idx_a(3,2)) ...
                == gid_b(idx_b(1,1):idx_b(1,2), idx_b(2,1):idx_b(2,2), idx_b(3,1):idx_b(3,2)));
        if sum(mask(:)) > most_match
            most_match = sum(mask(:));
            shift_best = shift_tmp;
        end
    end
    for j = 1:length(y_range)
        shift_tmp = shift_best;
        shift_tmp(2) = y_range(j);
        [idx_a, idx_b] = adjustSearchRange(shift_tmp, dim_a, dim_b);
        mask = (gid_a(idx_a(1,1):idx_a(1,2), idx_a(2,1):idx_a(2,2), idx_a(3,1):idx_a(3,2)) ...
                == gid_b(idx_b(1,1):idx_b(1,2), idx_b(2,1):idx_b(2,2), idx_b(3,1):idx_b(3,2)));
        if sum(mask(:)) > most_match
            most_match = sum(mask(:));
            shift_best = shift_tmp;
        end
    end
    for k = 1:length(z_range)
        shift_tmp = shift_best;
        shift_tmp(3) = z_range(k);
        [idx_a, idx_b] = adjustSearchRange(shift_tmp, dim_a, dim_b);
        mask = (gid_a(idx_a(1,1):idx_a(1,2), idx_a(2,1):idx_a(2,2), idx_a(3,1):idx_a(3,2)) ...
                == gid_b(idx_b(1,1):idx_b(1,2), idx_b(2,1):idx_b(2,2), idx_b(3,1):idx_b(3,2)));
        if sum(mask(:)) > most_match
            most_match = sum(mask(:));
            shift_best = shift_tmp;
        end
    end    
    
    % ----- Converge = the best shift didn't change -----
    if all(shift_best == shift_old)
        continue_opt = false;
    end
    
end

% disp('final')
% disp(shift)
end




