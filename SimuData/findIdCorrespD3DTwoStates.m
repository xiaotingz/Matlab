function corresp_d3d_45 = findIdCorrespD3DTwoStates(corresp_d3d_simu_an4, corresp_d3d_simu_an5, corresp_simu_new_45)
% ##################################################################
% * Notes
%     - IMPORTANT ASSUMPTION: there are only 1-to-1 and 1-to-2 corresps
%     in corresp_d3d_simu. This may not hold in general.
%     - Run in main_simuDataAnalysis.m
% * Input
%     - corresp_d3d_simu_an4 = {n1, 2}, corresp_d3d_simu_an5 = {n2, 2}
%           returned by findCorrespSimuD3D.m
%     - corresp_simu_new_45 = {n3, 2};
%           returned by findIdCorrespSimuTwoStates.m
% * Output
%     - corresp_d3d_45 = [id_d3d_an4, id_d3d_an5]
% ##################################################################
% -------------------------- differentiate the 1-to-1 and 1-to-2 correspondences --------------------------
grains_d3d_an4 = keys(corresp_d3d_simu_an4);
d3d_an4_multicorresp = [];
for key = 1:length(grains_d3d_an4)
    if length(corresp_d3d_simu_an4(key)) > 1
        d3d_an4_multicorresp = [d3d_an4_multicorresp; [key, corresp_d3d_simu_an4(key)]];
        remove(corresp_d3d_simu_an4, key);
    end
end
        
grains_d3d_an5 = keys(corresp_d3d_simu_an5);
d3d_an5_multicorresp = [];
for key = 1:length(grains_d3d_an5)
    if length(corresp_d3d_simu_an5(key)) > 1
        d3d_an5_multicorresp = [d3d_an5_multicorresp; [key, corresp_d3d_simu_an5(key)]];
        remove(corresp_d3d_simu_an5, key);
    end
end

% -------------------------- find corresp using the 1-to-1 correspondences --------------------------
tmp_an4 = [cell2mat(corresp_d3d_simu_an4.keys()); cell2mat(corresp_d3d_simu_an4.values())]';
tmp_an5 = [cell2mat(corresp_d3d_simu_an5.keys()); cell2mat(corresp_d3d_simu_an5.values())]';

corresp_d3d_45 = - ones(size(corresp_d3d_simu_an5, 1), 2);
id_d3d_an4 = cell2mat(corresp_d3d_simu_an4.keys());
id_simu_an4 = cell2mat(corresp_d3d_simu_an4.values());
id_simu_an5 = cell2mat(corresp_d3d_simu_an5.values());
corresp_simu_d3d_an5 = containers.Map(id_simu_an5, cell2mat(corresp_d3d_simu_an5.keys()));
id_simu_tracked_an4 = cell2mat(corresp_simu_new_45.keys())';

idx = 1;
for i = 1:length(id_d3d_an4)
    tmpid_d3d_an4 = id_d3d_an4(i);
    tmpid_simu_an4 = id_simu_an4(i);
    if ismember(tmpid_simu_an4, id_simu_tracked_an4)
        tmpid_simu_an5 = corresp_simu_new_45(tmpid_simu_an4);
        if ismember(tmpid_simu_an5, id_simu_an5)
            tmpid_d3d_an5 = corresp_simu_d3d_an5(tmpid_simu_an5);
            corresp_d3d_45(idx, :) = [tmpid_d3d_an4, tmpid_d3d_an5];
            idx = idx + 1;
        end
    end
end
corresp_d3d_45(idx:end, :) = [];

% -------------------------- find corresp using the 2-to-1 correspondences --------------------------
% """ 
% - NOTE in this case we only see 1-to-2 correspondences and it is assumed there 
% are only 1-to-1 and 2-to-1 corresps. This may not hold in general.
% - Outline
%     for each of id_d3d_an4: 
%         1) check its corresponding id_simu_new_an4s
%         2) find the corresponding id_simu_new_an5s
%         3) find the corresponding id_d3d_an5s
%     if only one unique id_d3d_an5 is found, build a mapping. Otherwise, error.
% """
% corresp_simu_new_45 = containers.Map(corresp_simu_new_45(:,1), corresp_simu_new_45(:,2));
corresp_simu_new_45_keys = cell2mat(keys(corresp_simu_new_45))';
for i = 1:size(d3d_an4_multicorresp)
    tmpid_d3d_an5_1 = -1;
    tmpid_d3d_an5_2 = -1;
    if ismember(d3d_an4_multicorresp(i, 2), corresp_simu_new_45_keys)
        tmpid_new_an4_1 = d3d_an4_multicorresp(i, 2);
        tmpid_new_an5_1= corresp_simu_new_45(tmpid_new_an4_1);
        [r_1, ~] = find(d3d_an5_multicorresp == tmpid_new_an5_1);
        tmpid_d3d_an5_1 = d3d_an5_multicorresp(r_1, 1);
    end
    if ismember(d3d_an4_multicorresp(i, 3), corresp_simu_new_45_keys)
        tmpid_new_an4_2 = d3d_an4_multicorresp(i, 3);
        tmpid_new_an5_2= corresp_simu_new_45(tmpid_new_an4_2);
        [r_2, ~] = find(d3d_an5_multicorresp == tmpid_new_an5_2);
        tmpid_d3d_an5_2 = d3d_an5_multicorresp(r_2, 1);
    end
    
    if tmpid_d3d_an5_1 > 0 && tmpid_d3d_an5_2 > 0
        if tmpid_d3d_an5_1 == tmpid_d3d_an5_2
            corresp_d3d_45 = [corresp_d3d_45; [d3d_an4_multicorresp(i, 1), tmpid_d3d_an5_1]];
        else
            warning(['MULTI corresps for corresp_d3d_45!']);
        end
    elseif tmpid_d3d_an5_1 > 0 || tmpid_d3d_an5_2 > 0
        corresp_d3d_45 = [corresp_d3d_45; [d3d_an4_multicorresp(i, 1), max(tmpid_d3d_an5_1, tmpid_d3d_an5_2)]];
    end
end

end