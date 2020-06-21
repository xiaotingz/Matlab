function [da_len_weighted, da_num_weighted] = calcGrainFaceDAs(faces, tl, tl_info)
% ############################################################################
% * Input
%   - faces = [n, 2]
%       The faces of interest. !!! The LEFT and RIGHT ORDER of the labels is important. !!! 
%   - tl = [m, 3]
%       Label of triple junctions, returned from findTripleLines.m
%   - tl_info = [m, 4]
%       [:, 1:3] = Turning angles at the Triple junctions. The angles are sorted in
%       accordance with grain_id in tl.
%       [:, 4] = length of the triple line.
% * Output
%   - da = [n, 3]
%       dihedral_angle = [da_l, da_r, da_oppo]
%       output different values, weighted by #TLs or length(TLs), because
%       it can be seen from paraview that the length(TLs) are not really good.
% * Note
%   - The left and right order of fl are arbitary to start with, but one
%   chosen needs to be consistent across all energy_grad features. See featurePrepEnergyGrad.m 
%   - Illustration (note r & l are not physical, just a convention)
%         face of interest = [r, l]
%                 |
%         r, da_r |  l, da_l
%                / \
%               /   \
%            da_opposite
% ############################################################################
% ----------------------- load debug data -----------------------
% % clear
% load('/Volumes/XIAOTING/Ni/working/190421_Hsmooth_energy_grad.mat','tracked_uniqueface_an4', 'tracked_uniqueface_an5');
% load('/Volumes/XIAOTING/Ni/working/190421_Hsmooth_TJs.mat');
% % % file = '/Users/xiaotingzhong/Desktop/Datas/synthetic/180502_CubicSingleEquiaxedOut.dream3d';
% file = '/Volumes/XIAOTING/Ni/An5new6_Hsmooth.dream3d';
% % % '/Volumes/XIAOTING/Ni/An4new6_fixOrigin2_smooth.dream3d';
% % % '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% % % '/Volumes/XIAOTING/Ni/An5new6_Hsmooth.dream3d'
% faces = tracked_uniqueface_an5;
% tl = triple_line_an5;
% tl_info = tl_info_an5;
% clear triple_line_an4 triple_line_an5 tl_info_an4 tl_info_an5
% ---------------------------------------------------------------
tl_da = tl_info(:, 1:3);
tl_len = tl_info(:,4);
da_len_weighted = zeros(size(faces, 1), 3);
da_num_weighted = zeros(size(faces, 1), 3);

for i = 1:size(faces, 1)
    % ##### Find Connecting TLs #####
    mask_connect = (sum(tl == faces(i, 1) | tl == faces(i, 2), 2) == 2);
    connect_tl = tl(mask_connect, :);
    connect_tl_da = tl_da(mask_connect, :);
    connect_tl_len = tl_len(mask_connect, :);
    
    % ##### Find DA With The RIGHT CORRESP #####
    mask_l = (connect_tl == faces(i, 1));
    mask_r = (connect_tl == faces(i, 2));
    mask_opp = ~(mask_l | mask_r);
    % """ weight the dihedral angles by triple line length """ 
    da_len_weighted(i, 1) = sum(connect_tl_da(mask_l) .* connect_tl_len) / sum(connect_tl_len);
    da_len_weighted(i, 2) = sum(connect_tl_da(mask_r) .* connect_tl_len) / sum(connect_tl_len);
    da_len_weighted(i, 3) = sum(connect_tl_da(mask_opp) .* connect_tl_len) / sum(connect_tl_len);
    % """ weight the dihedral angles by #triangles """ 
    da_num_weighted(i, 1) = sum(connect_tl_da(mask_l)) / size(mask_l, 1);
    da_num_weighted(i, 2) = sum(connect_tl_da(mask_r)) / size(mask_r, 1);
    da_num_weighted(i, 3) = sum(connect_tl_da(mask_opp)) / size(mask_opp, 1);
end

end 

% %% ######################################### Check ######################################### 
% % """
% % plot in paraview, see if reasonable
% % """
% face_id = h5read(file, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
% face_id(1,:) = [];
% 
% idx = randi(size(faces, 1));
% face_id = sort(face_id, 2);
% obj_label = sort(faces(idx, :));
% [idx_face, ~] = find(face_id(:,1) == obj_label(1) & face_id(:,2) == obj_label(2));
% 
% mask_connect = (sum(tl == faces(idx, 1) | tl == faces(idx, 2), 2) == 2);
% connect_tl = tl(mask_connect, :);
% connect_tl_da = tl_da(mask_connect, :);
% connect_tl_len = tl_len(mask_connect, :);
% 
% disp(['idx = ', num2str(idx)]);
% disp('-------------------------------------------------------------')
% disp(['face_label = ', num2str(faces(idx, :))]);
% disp(['face_id in an4:   ', num2str(idx_face)])
% disp('-----')
% disp(['da = ', num2str(tl_da(idx, :))])
% for i = 1:size(connect_tl, 1)
%     disp(['triple_line = ', num2str(connect_tl(i, :))]);
%     disp(['dihedral_angles = ', num2str(connect_tl_da(i, :))]);
%     disp(['length = ', num2str(connect_tl_len(i))]);
% end
% disp('-------------------------------------------------------------')
% 
% 
% 







