function id_corresp = findIdCorrespSameStateOldNew(id_given, id_new)
% ##################################################################
% * Objective
%     - This function is used to be used in simuDataAnalysis.
% * Input
%     - id_given = [n1, n2, n3]
%           a matrix of grain_id given by simulation, one grain can include
%           two isolated parts.
%     - id_new = [n1, n2, n3]
%           one island one id.
% * Output
%     - id_corresp = {num_given_grains, 1}, num_given_grains is the unique grain_ids in id_given
% ##################################################################
size_x = size(id_new, 1);
size_y = size(id_new, 2);
size_z = size(id_new, 3);
seen = - ones(length(unique(id_given)), 1);
idx = 1;

id_corresp = containers.Map('KeyType','double','ValueType','any');
for i = 1:size_x
    for j = 1:size_y
        for k = 1:size_z
            if ~ ismember(id_given(i, j, k), seen)
                id_corresp(id_given(i, j, k)) = id_new(i, j, k);
                seen(idx) = id_given(i, j, k);
                idx = idx + 1;
            else
                if ~ ismember(id_new(i, j, k), id_corresp(id_given(i, j, k))) 
                   id_corresp(id_given(i, j, k)) = [id_corresp(id_given(i, j, k)), id_new(i, j, k)];
                end
            end
        end
    end
end

end
