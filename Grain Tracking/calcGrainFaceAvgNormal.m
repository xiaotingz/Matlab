% function mask_ctwin = checkIfCoherentTwin(file, face_to_calc, thres_g, thres_n)
% ############################################################################
% * Input 
%     - face_to_calc = [n, 2], labels of the grain faces of interest.
%     - thres_norm, the threshold for angle between face_avg_norm and [111]
%  * Output
%     - mask_ctwin = [n, 1]
%  * Notes
%     - Mainly designed to check if twins are coherent or not. Use together
%     with getFaceRFvecs.m
% ############################################################################
% % ------------------ load data for debug --------------------
% file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
% face_to_calc = tracked_uniqueface_an4;
% rfvecs = rfvecs_an4;
thres_g = 0.06;
thres_n = 10;
% %-----------------------------------------------------------
% ##### Load and Clean Data #####
normal_ctwin = [1,1,1]';
normal_ctwin = normal_ctwin/sqrt(3);
rfvec_twin = normal_ctwin'*tand(60/2);
facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_normal = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
EA = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles'))';
EA(1,:) = [];
mask_inner = all(facelabel > 0, 2);
facelabel = facelabel(mask_inner, :);
tri_normal = tri_normal(mask_inner, :);

% 
% % ##### Get the Family of 111 direction #####
% O = CrysSym;
% ctwin_norm_family = zeros(24, 3);
% for i = 1:size(O, 3)
%     ctwin_norm_family(i, :) = O(:,:,i)*normal_ctwin;
% end
% ctwin_norm_family = unique(ctwin_norm_family, 'rows');
% n_family = size(ctwin_norm_family, 1);
% 

% ##### First Filter Out non-twin Misorientations #####
rfvec_twin = repmat(rfvec_twin, size(face_to_calc, 1), 1);
% [rfvecs]  = getFaceRFvecs(file, face_to_calc);
mask_ctwin = vecnorm(rfvecs - rfvec_twin, 2, 2) < thres_g;


% ##### Then Check Plane Normals #####
angle_diff_list = [];
for i = 1:size(face_to_calc, 1)
    if mask_ctwin(i)
        % ##### Get Grain Face Avg Normal #####
        % ----- Get triangles on the objective face -----
        mask_objface = (facelabel(:,1) == face_to_calc(i, 1) & facelabel(:,2) == face_to_calc(i, 2) |...
            facelabel(:,1) == face_to_calc(2) & facelabel(:,2) == face_to_calc(1));
        mask_reverse = facelabel(mask_objface,1) > facelabel(mask_objface,2);

        % ----- Get triangle normals, adjust order -----
        facetri_normal = tri_normal(mask_objface, :);
        facetri_normal(mask_reverse, :) = - facetri_normal(mask_reverse, :);

        % ##### Check if Angle Within Threshold #####
        % """
        % - Because the way D3D arrange normal is not clear, check {111}
        % instead of (111).
        % - To be sure matrix resizing is working correctly, make a set of
        % seudo facetri_normal using (111) direction and check results. 
        % """
        ctwin_norm_family = reshape(ctwin_norm_family, 1, []);
        normal_ctwin_mat = repmat(ctwin_norm_family, size(facetri_normal, 1), 1);
        normal_ctwin_mat = reshape(normal_ctwin_mat, [], 3);
        facetri_normal_mat = repmat(facetri_normal, n_family, 1);
        angle_diffs = abs(atand(vecnorm(cross(facetri_normal_mat, normal_ctwin_mat, 2), 2, 2)./dot(facetri_normal_mat, normal_ctwin_mat, 2)));
        angle_diffs = reshape(angle_diffs, [], n_family);
        angle_diffs = min(angle_diffs, [], 2);
        angle_diff_list = [angle_diff_list; angle_diffs];
       
        if sum(angle_diffs < thres_n) < length(angle_diffs)*0.8
            mask_ctwin(i) = false;            
        end
    end
end
    
    
histogram(angle_diff_list,'Normalization','probability')
set(gca, 'FontSize', 18);
xlabel('angle between triangle normals & \{111\}', 'FontSize', 18)
print('plane normal of twins','-dpng','-r300')










