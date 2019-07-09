function [search_range_a, search_range_b] = adjustSearchRange(shift, dim_a, dim_b)
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
%     - It works as keep volume_a stable and move volume_b. Whether
%     dim_a or dim_b larger are both fine.
% ##################################################################
search_range_a = [[1,1,1]; dim_a]';
search_range_b = [[1,1,1]; dim_b]';
% % 
% % if shift(1) < 0
% %     search_range_a(1,1) = 1;
% %     search_range_a(1,2) = min(shift(1) + dim_a(1), dim_b(1));
% %     search_range_b(1,1) = 1;
% %     search_range_b(1,2) = min(shift(1) + dim_a(1), dim_b(1));
% % elseif shift(1) + 1 + dim_a(1) > dim_b(1)
% %     search_range_a(1,1) = dim_a(1) + shift(1) - dim_b(1) + 1;
% %     search_range_a(1,2) = dim_a(1);
% %     search_range_b(1,1) = shift(1) + 1;
% %     search_range_b(1,2) = dim_b(1);
% % else
% %     search_range_a(1,1) = 1;
% %     search_range_a(1,2) = dim_a(1);
% %     search_range_b(1,1) = shift(1) + 1;
% %     search_range_b(1,2) = shift(1) + dim_a(1);
% % end
% %     
% % if shift(2) < 0
% %     search_range_a(2,1) = 1;
% %     search_range_a(2,2) = min(shift(2) + dim_a(2), dim_b(2));
% %     search_range_b(2,1) = 1;
% %     search_range_b(2,2) = min(shift(2) + dim_a(2), dim_b(2));
% % elseif shift(2) + 1 + dim_a(2) > dim_b(2)
% %     search_range_a(2,1) = dim_a(2) + shift(2) - dim_b(2) + 1;
% %     search_range_a(2,2) = dim_a(2);
% %     search_range_b(2,1) = shift(2) + 1;
% %     search_range_b(2,2) = dim_b(2);
% % else
% %     search_range_a(2,1) = 1;
% %     search_range_a(2,2) = dim_a(2);
% %     search_range_b(2,1) = shift(2) + 1;
% %     search_range_b(2,2) = shift(2) + dim_a(2);
% % end
% % 
% % if shift(3) < 0
% %     search_range_a(3,1) = 1;
% %     search_range_a(3,2) = min(shift(3) + dim_a(3), dim_b(3));
% %     search_range_b(3,1) = 1;
% %     search_range_b(3,2) = min(shift(3) + dim_a(3), dim_b(3));
% % elseif shift(3) + 1 + dim_a(3) > dim_b(3)
% %     search_range_a(3,1) = dim_a(3) + shift(3) - dim_b(3) + 1;
% %     search_range_a(3,2) = dim_a(3);
% %     search_range_b(3,1) = shift(3) + 1;
% %     search_range_b(3,2) = dim_b(3);
% % else
% %     search_range_a(3,1) = 1;
% %     search_range_a(3,2) = dim_a(3);
% %     search_range_b(3,1) = shift(3) + 1;
% %     search_range_b(3,2) = shift(3) + dim_a(3);
% % end

if shift(1) < 0
    search_range_a(1,1) = 1;
    search_range_a(1,2) = min(dim_a(1), dim_b(1) + shift(1));
    search_range_b(1,1) = 1 - shift(1);
    search_range_b(1,2) = min(dim_a(1) - shift(1), dim_b(1));
else
    search_range_a(1,1) = 1 + shift(1);
    search_range_a(1,2) = min(dim_a(1), dim_b(1) + shift(1));
    search_range_b(1,1) = 1;
    search_range_b(1,2) = min(dim_a(1) - shift(1), dim_b(1));
end

if shift(2) < 0
    search_range_a(2,1) = 1;
    search_range_a(2,2) = min(dim_a(2), dim_b(2) + shift(2));
    search_range_b(2,1) = 1 - shift(2);
    search_range_b(2,2) = min(dim_a(2) - shift(2), dim_b(2));
else
    search_range_a(2,1) = 1 + shift(2);
    search_range_a(2,2) = min(dim_a(2), dim_b(2) + shift(2));
    search_range_b(2,1) = 1;
    search_range_b(2,2) = min(dim_a(2) - shift(2), dim_b(2));
end

if shift(3) < 0
    search_range_a(3,1) = 1;
    search_range_a(3,2) = min(dim_a(3), dim_b(3) + shift(3));
    search_range_b(3,1) = 1 - shift(3);
    search_range_b(3,2) = min(dim_a(3) - shift(3), dim_b(3));
else
    search_range_a(3,1) = 1 + shift(3);
    search_range_a(3,2) = min(dim_a(3), dim_b(3) + shift(3));
    search_range_b(3,1) = 1;
    search_range_b(3,2) = min(dim_a(3) - shift(3), dim_b(3));
end

end
    
    
% ####################################### Test #######################################
% file_a = file_an3;
% file_b = file_an4;
% dim_a = h5read(file_a, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS');
% dim_b = h5read(file_b, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS');
% 
% % """
% % Initial search range:
% %       Note volume_a stable and volume_b moving, 
% %       If dim_a > dim_b, should do more positive motion. 
% %           volume_a    |--------------------------|
% %           volume_b  - |----------------| +                        (mobile)
% %       If dim_a < dim_b, should do more negative motion.
% %           volume_a    |--------------------------|
% %           volume_b  - |------------------------------------| +    (mobile)
% % """
% ranges = cell(3, 1);
% for i = 1:3
%     diff = abs(dim_b(i) - dim_a(i));
%     if dim_a(i) > dim_b(i)
%         ranges{i} = (-3 : diff + 3)';
%     else
%         ranges{i} = (-3 - diff :  3)';
%     end
% end
% 
% 
% disp('###################################################')
% disp([ranges{1}(1), ranges{1}(end); ranges{2}(1), ranges{2}(end); ranges{3}(1), ranges{3}(end)])
% disp([0,0,0])
% [search_range_a, search_range_b] = adjustSearchRange([0,0,0], dim_a', dim_b');
% disp([search_range_a, [-1; -1; -1], search_range_b])
% disp([-1,-1,-1])
% [search_range_a, search_range_b] = adjustSearchRange([-1,-1,-1], dim_a', dim_b');
% disp([search_range_a, [-1; -1; -1], search_range_b])
% disp([1,1,1])
% [search_range_a, search_range_b] = adjustSearchRange([1,1,1], dim_a', dim_b');
% disp([search_range_a, [-1; -1; -1], search_range_b])
% disp([4, -1, 3])
% [search_range_a, search_range_b] = adjustSearchRange([4, -1, 3], dim_a', dim_b');
% disp([search_range_a, [-1; -1; -1], search_range_b])
% disp('###################################################')
% 
% for i = 1:10
%     rng('shuffle')
%     shift(1) = ranges{1}(randi(length(ranges{1})));
%     shift(2) = ranges{2}(randi(length(ranges{2})));
%     shift(3) = ranges{3}(randi(length(ranges{3})));
%     [search_range_a, search_range_b] = adjustSearchRange(shift, dim_a', dim_b');
%     volume_wrong = (search_range_b(1,2) - search_range_b(1,1)) ~= (search_range_a(1,2) - search_range_a(1,1)) || ...
%         (search_range_b(2,2) - search_range_b(2,1)) ~= (search_range_a(2,2) - search_range_a(2,1)) || ...
%         (search_range_b(3,2) - search_range_b(3,1)) ~= (search_range_a(3,2) - search_range_a(3,1));
%     idx_bad = search_range_b(1,2) > dim_b(1) || search_range_a(1,2) > dim_a(1) || ...
%         search_range_b(2,2) > dim_b(2) || search_range_a(2,2) > dim_a(2) || ...
%         search_range_b(3,2) > dim_b(3) || search_range_a(3,2) > dim_a(3);
%     if volume_wrong || idx_bad
%         if volume_wrong
%             warning(['volume_wrong at: ', num2str(shift)]);
%         else
%             warning(['idx_bad at: ', num2str(shift)]);
%         end
%     end
%     disp('-----------------')
%     disp(shift)
%     disp([search_range_a, [-1; -1; -1], search_range_b])
% end
% 
% 



    