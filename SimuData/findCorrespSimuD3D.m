function corresp_d3d_simu = findCorrespSimuD3D(file, id_new)
% ##################################################################
% * Notes
%     - Finds corresp between grain_id_cleaned_simulation and grain_id_d3d,
%     for the same state.
%     - Dependency: findIdCorrespSameStateOldNew.m
% * Input
%     - file
%           D3D reconstructed file.
%     - id_new_ = [n1, n2, n3]
%           returned by findNewId.m and cleaned by findSizesAndCentroids.m
% * Output
%     - corresp_d3d_simu = {grain_id_d3d, grain_id_new}
% ##################################################################
% % % ----------------------- load debug data -----------------------
% load('simu_data.mat', 'id_newclean_an5');
% id_new = id_newclean_an5;
% clear id_newclean_an5
% file = '/Volumes/XIAOTING/Ni/simu_An5_clean_seg.dream3d';
% % % ---------------------------------------------------------------

% ----------------------- read data -----------------------
id_d3d = squeeze(h5read(file, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));
id_new = int32(id_new);
disp([length(unique(id_d3d)), length(unique(id_new))]);

% ----------------------- read data -----------------------
corresp_d3d_simu = findIdCorrespSameStateOldNew(id_d3d, id_new);
remove(corresp_d3d_simu, 0);


keys_multi_corresp = [];
all_keys = keys(corresp_d3d_simu)';
for i = 1:length(all_keys)
    if ismember(-1, corresp_d3d_simu(all_keys{i}))
        corresp_d3d_simu(all_keys{i}) = setdiff(corresp_d3d_simu(all_keys{i}), -1);
    end
    if length(corresp_d3d_simu(all_keys{i})) > 1
        keys_multi_corresp = [keys_multi_corresp; all_keys{i}];
    end
end

end












