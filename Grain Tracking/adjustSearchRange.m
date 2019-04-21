function [search_range_an4, search_range_an5] = adjustSearchRange(shift, dim_an4, dim_an5)
% ##################################################################
% * Input
%     - shift = [3, 1] = [shift_x, shift_y, shift_z]
%     - dim_ = [3, 1], dimension of the volume
% * Output
%     - search_range = [3, 2]
%         [x_start, x_end; y_start, y_end; z_start, z_end]
% * Note
%     - This script is to be used in calcVolumeDisplacement.m 
%     - The purpose is to keep the two search volumes legal after applying a shift.
% ##################################################################
search_range_an4 = [[1,1,1]; dim_an4]';
search_range_an5 = [[1,1,1]; dim_an4]';

if shift(1) < 0
    search_range_an4(1,1) = 1;
    search_range_an4(1,2) = shift(1) + dim_an4(1);
    search_range_an5(1,1) = 1;
    search_range_an5(1,2) = shift(1) + dim_an4(1);
elseif shift(1) + 1 + dim_an4(1) > dim_an5(1)
    search_range_an4(1,1) = dim_an4(1) + shift(1) - dim_an5(1) + 1;
    search_range_an4(1,2) = dim_an4(1);
    search_range_an5(1,1) = shift(1) + 1;
    search_range_an5(1,2) = dim_an5(1);
else
    search_range_an4(1,1) = 1;
    search_range_an4(1,2) = dim_an4(1);
    search_range_an5(1,1) = shift(1) + 1;
    search_range_an5(1,2) = shift(1) + dim_an4(1);
end
    
if shift(2) < 0
    search_range_an4(2,1) = 1;
    search_range_an4(2,2) = shift(2) + dim_an4(2);
    search_range_an5(2,1) = 1;
    search_range_an5(2,2) = shift(2) + dim_an4(2);
elseif shift(2) + 1 + dim_an4(2) > dim_an5(2)
    search_range_an4(2,1) = dim_an4(2) + shift(2) - dim_an5(2) + 1;
    search_range_an4(2,2) = dim_an4(2);
    search_range_an5(2,1) = shift(2) + 1;
    search_range_an5(2,2) = dim_an5(2);
else
    search_range_an4(2,1) = 1;
    search_range_an4(2,2) = dim_an4(2);
    search_range_an5(2,1) = shift(2) + 1;
    search_range_an5(2,2) = shift(2) + dim_an4(2);
end

if shift(3) < 0
    search_range_an4(3,1) = 1;
    search_range_an4(3,2) = shift(3) + dim_an4(3);
    search_range_an5(3,1) = 1;
    search_range_an5(3,2) = shift(3) + dim_an4(3);
elseif shift(3) + 1 + dim_an4(3) > dim_an5(3)
    search_range_an4(3,1) = dim_an4(3) + shift(3) - dim_an5(3) + 1;
    search_range_an4(3,2) = dim_an4(3);
    search_range_an5(3,1) = shift(3) + 1;
    search_range_an5(3,2) = dim_an5(3);
else
    search_range_an4(3,1) = 1;
    search_range_an4(3,2) = dim_an4(3);
    search_range_an5(3,1) = shift(3) + 1;
    search_range_an5(3,2) = shift(3) + dim_an4(3);
end

end
    
    
% ####################################### Test #######################################
% [tmp_an4, tmp_an5] = adjustSearchRange([6,0,0], dim_an4, dim_an5);
% disp(['tmp_an4_x = [', num2str(tmp_an4(1,:)), ']', ',   diff = ', num2str(tmp_an4(1,2) - tmp_an4(1,1))]);
% disp(['tmp_an5_x = [', num2str(tmp_an5(1,:)), ']', ',   diff = ', num2str(tmp_an5(1,2) - tmp_an5(1,1))]);
    