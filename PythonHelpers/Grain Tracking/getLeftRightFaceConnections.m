% function nn_faces = getLeftRightFaceConnections(unique_faces, tls, file)
% ##########################################################################
% * Input
%     - unique_faces = [n1,2]  
%           returned by trackFace.m or trackUniqueFace.m, usually named tracked_uniqueface_ in other functions
%     - file: file name
%     - tls = [m, 3]
%           returned by Topologies/findTripleLines.m
% * Output
%     - nn_faces = {{k1, 2}*n, {k2, 2}*n}
%           k1 is #left_connected_faces, k2 is #right_connected_faces, of the ith grain face. 
%           There are in total n grain faces.
% * NOTE
%     - This function is to calculate left/right connection of grain 
%       faces. The idea is related to that of high/low energy (big/small size) grains, 
%       but applied on the level of grain faces. The objective is to predict migration sign.
%     - Related functions: calcFaceToGrainCentroidDist.m
% ##########################################################################
% % ----------------------- load debug data -----------------------
tracked_faces = tracked_uniqueface_an5;
tls = tl_an5;
% clear tracked_uniqueface_an4 tl_an4
file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
% % ---------------------------------------------------------------
all_faces = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
all_faces = unique(all_faces(all(all_faces>0, 2), :), 'rows');
all_faces = sort(all_faces, 2);

nn_faces = cell(size(tracked_faces, 1), 2);
for i = 1:size(tracked_faces, 1)
    left_grain = tracked_faces(i, 1);
    right_grain = tracked_faces(i, 2);
    mask_tls = (sum((tls == left_grain | tls == right_grain), 2) == 2);
    tl = tls(mask_tls, :);
    candidate_grains = setdiff(unique(tl), [left_grain, right_grain]);
    for j = 1:size(tl, 1)
        candidate = candidate_grains(j);
        if ismember(sort([candidate, left_grain]), all_faces)
            nn_faces{i, 1} = [nn_faces{i, 1}; sort([candidate, left_grain])];
        end
        if ismember(sort([candidate, right_grain]), all_faces)
            nn_faces{i, 2} = [nn_faces{i, 2}; sort([candidate, right_grain])];
        end
    end
end









