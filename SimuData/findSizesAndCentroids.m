function [sizes, centroids, id] = findSizesAndCentroids(id)
% ##################################################################
% * Objective
%     - This function is used to be used in simuDataAnalysis.
% * Input
%     - id = [n1, n2, n3], a matrix of grain_id
%          given by findNewId.m
% * Output
%     - sizes = [num_grains, 1], num_grains is the unique grain_ids in id
%     - centroids = [num_grains, 3]
% ##################################################################
% % ----------------------- load debug data -----------------------
% id = id_new_an4;
% % ---------------------------------------------------------------

size_x = size(id, 1);
size_y = size(id, 2);
size_z = size(id, 3);

% ----- cal position indexes -----
[r, c] = find(ones(size_x, size_y));
x_idxes = repmat(reshape(r, [size_x, size_y]), 1, 1, size_z);
y_idxes = repmat(reshape(c, [size_x, size_y]), 1, 1, size_z);
z_idxes = ones(size(id));
for k = 1:size(z_idxes, 3)
    z_idxes(:,:,k) = z_idxes(:,:,k) * k;
end

% ----- position indexes -----
id_unique = unique(id);
sizes = zeros(size(id_unique));
centroids = zeros(size(id_unique, 1), 3);
thres = 27;

for i = 1 : length(id_unique)
    mask = (id == id_unique(i));
    
    sizes(i) = sum(mask(:));
    if sizes(i) < thres
        id(mask) = -1;
    end
    
    x = x_idxes(mask);
    y = y_idxes(mask);
    z = z_idxes(mask);
    centroids(i, :) = [sum(x(:)), sum(y(:)), sum(z(:))] ./ sizes(i);
end

end