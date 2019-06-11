function id_new = findNewId(id_given)
% ##################################################################
% * Objective
%     - Some grains returned by the simulation contain two different
%     islands. This scipt reassigns grain_ids s.t. every island has a
%     unique grain_id
% * Input
%     - id_given = [n1, n2, n3], a matrix of grain_id given by simulation
% * Output
%     - id_new = [n1, n2, n3], a matrix of reassigned grain_id
% ##################################################################
% % ----------------------- load debug data -----------------------
% id_new = ones(3, 3, 3);
% id_new(1,1,1) = 2;
% id_new(2,2,2) = 3;
% % ---------------------------------------------------------------

id_new = id_given;
new = ones(size(id_new));
size_x = size(id_new, 1);
size_y = size(id_new, 2);
size_z = size(id_new, 3);

cnt = 0 ;
for i = 1:size_x
    for j = 1:size_y
        for k = 1:size_z
            if new(i, j, k) > 0
                cnt = cnt + 1;
                [id_new, new] = exploreIsland(id_new, new, i, j, k, cnt);
                disp(cnt);
            end
        end
    end
end

end