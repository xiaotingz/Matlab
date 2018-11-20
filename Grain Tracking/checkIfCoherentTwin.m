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
file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
face_to_calc = tracked_uniqueface_an4;
rfvecs = rfvecs_an4;
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
EA = rad2deg(EA);
mask_inner = all(facelabel > 0, 2);
facelabel = facelabel(mask_inner, :);
tri_normal = tri_normal(mask_inner, :);


% ##### First Filter Out non-twin Misorientations #####
rfvec_twin = repmat(rfvec_twin, size(face_to_calc, 1), 1);
% [rfvecs]  = getFaceRFvecs(file, face_to_calc);
mask_ctwin = vecnorm(rfvecs - rfvec_twin, 2, 2) < thres_g;
%%

% ##### Then Check Plane Normals #####
angle_diff_list = [];
angle_diff_cnt = [];
O = CrysSym;
% """
% Note the symmetry matrixes are transposed this way. But the group should
% stay the same. 
% """
for i = 1:size(face_to_calc, 1)
    disp(i)
    if mask_ctwin(i)
        % ----- Get triangle normal in sample frame -----
        mask_objface = (facelabel(:,1) == face_to_calc(i, 1) & facelabel(:,2) == face_to_calc(i, 2) |...
            facelabel(:,1) == face_to_calc(2) & facelabel(:,2) == face_to_calc(1));
        facetri_normal = tri_normal(mask_objface, :);
        angle_diffs = ones(size(facetri_normal, 1), 1)*90;

        % ----- Convert to crystal frame and apply symmetries -----
        g1 = EAtoG(EA(face_to_calc(i, 1), :));
        g2 = EAtoG(EA(face_to_calc(i, 1), :));
        for j = 1:size(facetri_normal, 1)
            normal = facetri_normal(j, :)';
            for k = 1:size(O, 3)
                for l = 1:size(O, 3)
                    normal_vari = (O(:,:,k) * g1) * normal;
                    normal_ctwin_vari = (O(:,:,l) * g1) * normal_ctwin;
                    ang_tmp = abs(atand(norm(cross(normal_vari, normal_ctwin_vari))/dot(normal_vari, normal_ctwin_vari)));
                    if ang_tmp < angle_diffs(j)
                        angle_diffs(j) = ang_tmp;
                    end
                    normal_vari = O(:,:,k) * g2 * normal;
                    ang_tmp = abs(atand(norm(cross(normal_vari, normal_ctwin_vari))/dot(normal_vari, normal_ctwin_vari)));
                    if ang_tmp < angle_diffs(j)
                        angle_diffs(j) = ang_tmp;
                    end
                end
            end
        end
        
        angle_diff_list = [angle_diff_list; angle_diffs];
        angle_diff_cnt = [angle_diff_cnt; sum(mask_objface)];
        
        % ##### Check if Angle Within Threshold #####
        if sum(angle_diffs < thres_n) < length(angle_diffs)*0.8
            mask_ctwin(i) = false;            
        end
    end
end
    
    
histogram(angle_diff_list,'Normalization','probability')
set(gca, 'FontSize', 18);
xlabel('angle between triangle normals & \{111\}', 'FontSize', 18)
print('plane normal of twins','-dpng','-r300')









