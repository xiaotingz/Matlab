function neighbors = getNeighList(grain_id, num_neigh, neighbor_list)
% ############################################################################
% NOTES
% - assuming the non-function 0 at the first element of numNeighbors has been deleted.
% ############################################################################
if num_neigh(1) == 0
    warning('Delete the first element in num_neigh before running this script!')
end
neighbors =  neighbor_list(sum(num_neigh(1:grain_id-1)) + 1 : sum(num_neigh(1:grain_id)));
end