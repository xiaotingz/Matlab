function [mat, new] = exploreIsland(mat, new, i, j, k, new_id)
% ##################################################################
% * Objective
%     - This function is used to correcting simulation grain ID. Some
%     simulation grains have two isolated pieces with the same grain ID.
%     We'll assign grains to a different set of IDs in which every isolated
%     piece have a different ID. 
% * Input
%     - mat = [n1, n2, n3], a matrix of grain_id given by simulaiton. 
%         every time exploreIsland is called, all voxels of a single grain
%         will be visited, and the grain_id will be changed to new_id.
%     - new = [n1, n2, n3], indicating if this position has been visited.
%     - i, j, k = scalar, position indexes.
% ##################################################################
% % ----------------------- load debug data -----------------------
% mat = ones([2,2,3]);
% new = ones(size(mat));
% i = 1;
% j = 1;
% k = 1;
% mat(2, 2, 3) = 2;
% % ---------------------------------------------------------------

next_pos = [i, j, k];
grain_id = mat(i, j, k);

while ~ isempty(next_pos)
    i = next_pos(1, 1);
    j = next_pos(1, 2);
    k = next_pos(1, 3);
    next_pos(1, :) = [];
    legal_pos = (0 < i && i <= size(mat, 1) && 0 < j && j <= size(mat, 2) && 0 < k && k <= size(mat, 3));
    explore = legal_pos && new(i, j, k) > 0 && mat(i, j, k) == grain_id;
    if explore
        new(i, j, k) = 0;
        mat(i, j, k) = new_id;
        idx = size(next_pos, 1);
%         next_pos(idx+1, :) = [i+1, j, k];
%         next_pos(idx+2, :) = [i, j+1, k];
%         next_pos(idx+3, :) = [i, j, k+1];
%         next_pos(idx+4, :) = [i-1, j, k];
%         next_pos(idx+5, :) = [i, j-1, k];
%         next_pos(idx+6, :) = [i, j, k-1];
        next_pos(idx+1 : idx+6, :) = [[i+1, j, k]; [i, j+1, k]; [i, j, k+1]; ...
                                      [i-1, j, k]; [i, j-1, k]; [i, j, k-1]];
    end
end

end

