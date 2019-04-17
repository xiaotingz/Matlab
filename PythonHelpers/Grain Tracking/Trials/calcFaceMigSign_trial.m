% function mig_signs = calcFaceMigSign(file_an4, file_an5, tracked_uniqueface_an4, tracked_uniqueface_an5, grain_data_an4, grain_data_an5)
% ##########################################################################
% * Input
%     - tracked_uniqueface_ = [n,2]  
%           returned by trackFace.m or trackUniqueFace.m
%     - grain_data_ = [m, k]
%           m = #grains, k = #energies passed in
%           size, F, or F - <Fnn>, or integral curvature, grain area
% * Output
%     - mig_signs = [n, k]
%           n = #grain_faces, k = #energies passed in
% * NOTE
%     - This function is to calculate the sign of face migration. 
%     - sign = 1    if face moving towards the grain with lower grain_data; 
%       sign = -1   if face moving toward the grain with higher grain_data;
%       sign = 0    if the two grains swapped grain_data level.
%     - Dependence: calcFaceToGrainCentroidDist.m

% """
% TODO: there are two ways to define grain_data relavance: just an4 only, or
% use grain_data difference. now is using an4 only. 
% """
% ##########################################################################
% ----------------------- load debug data -----------------------
% file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
% file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/181107.mat', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5');
% load('look_up_table_an4_an5.mat')

% ##### run  /Grain Curvature/G_F_mF.m to get the following data: from data_grain & F_mF_diff #####
% grain_data_an4 = {grain_size_an4, grain_F_an4, grain_itg_curv_an4, F_avgFnn_an4};
% grain_data_an5 = {grain_size_an5, grain_F_an5, grain_itg_curv_an5, F_avgFnn_an5};
% ---------------------------------------------------------------


%% #################################### Caculate sign, grain properties ####################################
% ##### Recalc tracked_uniqueface_an4 To Keep Grains Consistent With track_uniquefaces_an5 #####
% ----- make a real index table out of look_up_table -----
look_up_table = sortrows(look_up_table, 2);
look_up_table_2to1 = zeros(max(look_up_table(:,2)),2);
idx = 1;
for i = 1 : max(look_up_table(:,2))
    if look_up_table(idx,2) == i
        look_up_table_2to1(i,1) = look_up_table(idx,1);
        look_up_table_2to1(i,2) = look_up_table(idx,2);
        idx = idx + 1;
    else
        look_up_table_2to1(i,1) = NaN;
        look_up_table_2to1(i,2) = i;
    end
end
% ----- remake tracked_uniqueface_an4 -----
tracked_uniqueface_an4 = zeros(size(tracked_uniqueface_an5));
for i = 1:size(tracked_uniqueface_an5, 1)
    for j = 1:size(tracked_uniqueface_an5, 2)
        tracked_uniqueface_an4(i, j) = look_up_table_2to1(tracked_uniqueface_an5(i,j), 1);
    end
end

% #################################### Calc Face To Grain Centroid Dist ####################################
% ----- Absolute distance -----
% dists_f_g_an4 = calcFaceToGrainCentroidDist(file_an4, tracked_uniqueface_an4);
% dists_f_g_an5 = calcFaceToGrainCentroidDist(file_an5, tracked_uniqueface_an5);
% ----- Convert distance to portion -----
total_dists_an4 = sum(dists_f_g_an4, 2);
total_dists_an5 = sum(dists_f_g_an5, 2);
dists_f_g_an4 = dists_f_g_an4 ./ total_dists_an4;
dists_f_g_an5 = dists_f_g_an5 ./ total_dists_an5;

% % ##### Check Sign For Each Grain grain_data Metric #####
eps = 1e-2;
sign = zeros(size(tracked_uniqueface_an4, 1), size(grain_data_an4, 2));
for i = 1:size(grain_data_an4, 2)
    for j = 1:size(tracked_uniqueface_an4, 1)
        grain_data_an4_cur = grain_data_an4{i};
        grain_data_an5_cur = grain_data_an5{i};
        id_l_an4 = tracked_uniqueface_an4(j, 1);
        id_r_an4 = tracked_uniqueface_an4(j, 2);
        id_l_an5 = tracked_uniqueface_an5(j, 1);
        id_r_an5 = tracked_uniqueface_an5(j, 2);
        % """
        % eps: a small variable to filter out numerical errors.
        % smaller size = bigger grain_data ???????????
        % only grain_data in an4 should be considered, or grain_data in an4 & an5 both
        % needs to be considered ???????????
        % """
%         left_grain_grain_data_big = (grain_data_an4_cur(id_l_an4) < grain_data_an4_cur(id_r_an4) && ...
%                                 grain_data_an5_cur(id_l_an5) < grain_data_an5_cur(id_r_an5));
    %     right_grain_grain_data_big = (grain_data_an4_cur(id_l_an4) > grain_data_an4_cur(id_r_an4) && ...
    %                             grain_data_an5_cur(id_l_an5) > grain_data_an5_cur(id_r_an5));
        left_grain_small = (grain_data_an4_cur(id_l_an4) < grain_data_an4_cur(id_r_an4));
        if left_grain_small
            if dists_f_g_an4(j, 1) < (dists_f_g_an5(j, 1) - eps)
                sign(j, i) = -1;
            elseif dists_f_g_an4(j, 1) > (dists_f_g_an5(j, 1) + eps)
                sign(j, i) = 1;
            else
                sign(j, i) = 0;
            end
        else
            if dists_f_g_an4(j, 1) < (dists_f_g_an5(j, 1) - eps)
                sign(j, i) = 1;
            elseif dists_f_g_an4(j, 1) > (dists_f_g_an5(j, 1) + eps)
                sign(j, i) = -1;
            else
                sign(j, i) = 0;
            end
        end

    end
end

% ----- Adjust sign for itg_grain_curv -----
% """
% small size = small #faces = large itg_grain_curv = small F-<Fnn>
% """
sign(:,3) = -sign(:,3);

%%
% % #################################### Plot Sign v.s. Misorientation Angle ####################################
% ##### Calculate Misorientation Angles #####
mis_ang_full = double(h5read(file_an4,'/DataContainers/ImageDataContainer/CellFeatureData/MisorientationList'))';
face_full_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';

sign = sign_areadiff;

mask = all(face_full_an4>0, 2);
face_full_an4 = face_full_an4(mask, :);
mis_ang_an4 = zeros(size(sign,1), 1);
for i = 1:size(tracked_uniqueface_an4, 1)
    mask = (face_full_an4(:,1) == tracked_uniqueface_an4(i,1) & face_full_an4(:,2) == tracked_uniqueface_an4(i,2) | ...
            face_full_an4(:,1) == tracked_uniqueface_an4(i,2) & face_full_an4(:,2) == tracked_uniqueface_an4(i,1));
    mis_ang_an4(i) = sum(mis_ang_full(mask))/sum(mask);
end


% ##### Bin Misorientation Angles And Signs #####
bin_edges = (0:10:60)';
mis_ang_an4_binned = zeros(size(bin_edges,1),1);
sign_total_binned = zeros(size(bin_edges,1), size(sign, 2));
sign_pos_binned = zeros(size(bin_edges,1), size(sign, 2));
sign_neg_binned = zeros(size(bin_edges,1), size(sign, 2));
sign_zero_binned = zeros(size(bin_edges,1), size(sign, 2));
for i = 1:size(mis_ang_an4,1)
    idx = floor(mis_ang_an4(i)/10) + 1;
    mis_ang_an4_binned(idx) = mis_ang_an4_binned(idx) + 1;
    for k = 1:size(sign, 2)
        sign_total_binned(idx, k) = sign_total_binned(idx, k) + sign(idx, k);
        if sign(i, k) > 0
            sign_pos_binned(idx, k) = sign_pos_binned(idx, k) + 1;
        elseif sign(i, k) < 0
            sign_neg_binned(idx, k) = sign_neg_binned(idx, k) + 1;
        else
            sign_zero_binned(idx, k) = sign_zero_binned(idx, k) + 1;
        end
    end
end
      
%%
set(0,'defaultAxesLabelFontSize',1.2)
set(0,'defaultAxesFontSize',18)
figure
% scatter(bin_edges, sign_pos_binned(:,1))
idx_grain_data = 2;
bar(bin_edges, [sign_pos_binned(:,idx_grain_data), sign_neg_binned(:,idx_grain_data), sign_zero_binned(:,idx_grain_data)])
xlabel('misorientation angle, ^{\circ}')
ylabel('#Faces')
legend('pos motion', 'neg motion', 'no motion', 'Location', 'northwest')
fig_name = ('grain_data gradient from max(\Delta NNF\_area), eps = ' + string(eps));
title(fig_name)
print(fig_name,'-dpng','-r300')





%% #################################### Calculate Sign, Face Properties ####################################
% % ##### Calculate face areas #####
% """ faceCurvs = [N,2] = [integralArea, integralCurvature] """
% [faces_an4, faces_an5, face_corresp] = trackFace(file_an4, file_an5, look_up_table, 'use_complete_faces');
% face_itg_curv_an4_all = calcFaceItgCurv(file_an4, faces_an4, 'all_faces');
% face_itg_curv_an5_all = calcFaceItgCurv(file_an5, faces_an5, 'all_faces');
% mask_an4 = (faces_an4(:,1) < faces_an4(:,2));
% mask_an5 = (faces_an5(:,1) < faces_an5(:,2));
% faces_an4_all = faces_an4(mask_an4, :);
% faces_an5_all = faces_an5(mask_an5, :);
% face_itg_curv_an4_all = face_itg_curv_an4_all(mask_an4, :);
% face_itg_curv_an5_all = face_itg_curv_an5_all(mask_an5, :);
% tracked_uniqueface_an4_ordered = sort(tracked_uniqueface_an4, 2);


% % ##### Calculate face areas change for all faces in an4 #####
% """
%   - Note that in this logic, only face area decrease has been considered.
%     face area increase doesn't matter. 
%   - May need to consider face area change due to field of view change. 
% """
area_diff_an4 = - face_itg_curv_an4_all(:, 1);
mask_tracked_an4 = ismember(faces_an4_all, tracked_uniqueface_an4_ordered, 'rows');
for i = 1:size(faces_an4_all, 1)
    if mask_tracked_an4(i) == 1
        mask = ismember(tracked_uniqueface_an4_ordered, faces_an4_all(i, :), 'rows');
        f_area_an5 = face_itg_curv_an5(mask, 1);
        % --- area_diff_an4 has been initialized to be negative ---
        area_diff_an4(i) = f_area_an5 + area_diff_an4(i);
    end
end
        
%% ############################## Energy gradient from neighboring faces ##############################
sign_areadiff = zeros(size(tracked_uniqueface_an4, 1), 2);
fnnf_maxadec_l_m_r = zeros(size(tracked_uniqueface_an4, 1), 1);
for i = 1:length(sign_areadiff)
    left_nnf_an4 = nn_faces_an4{i, 1};
    mask_lnnf_an4 = ismember(faces_an4_all, left_nnf_an4, 'rows');
    left_nnf_area_maxdec_an4 = min(area_diff_an4(mask_lnnf_an4, :));
    left_nnf_area_avgdec_an4 = sum(area_diff_an4(mask_lnnf_an4, :))/sum(mask_lnnf_an4);
    right_nnf_an4 = nn_faces_an4{i, 2};
    mask_rnnf_an4 = ismember(faces_an4_all, right_nnf_an4, 'rows');
    right_nnf_area_maxdec_an4 = min(area_diff_an4(mask_rnnf_an4, :));
    right_nnf_area_avgdec_an4 = sum(area_diff_an4(mask_rnnf_an4, :))/sum(mask_rnnf_an4);
    
    fnnf_maxadec_l_m_r(i) = left_nnf_area_maxdec_an4 - right_nnf_area_maxdec_an4;
 
    left_nnf_shrink_more = (left_nnf_area_maxdec_an4 < right_nnf_area_maxdec_an4);
    if left_nnf_shrink_more
        if dists_f_g_an4(i, 1) < (dists_f_g_an5(i, 1) - eps)
            sign_areadiff(i, 1) = -1;
        elseif dists_f_g_an4(i, 1) > (dists_f_g_an5(i, 1) + eps)
            sign_areadiff(i, 1) = 1;
        else
            sign_areadiff(i, 1) = 0;
        end
    else
        if dists_f_g_an4(i, 1) < (dists_f_g_an5(i, 1) - eps)
            sign_areadiff(i, 1) = 1;
        elseif dists_f_g_an4(i, 1) > (dists_f_g_an5(i, 1) + eps)
            sign_areadiff(i, 1) = -1;
        else
            sign_areadiff(i, 1) = 0;
        end
    end
    
    left_nnf_shrink_more = (left_nnf_area_avgdec_an4 < right_nnf_area_avgdec_an4);
    if left_nnf_shrink_more
        if dists_f_g_an4(i, 1) < (dists_f_g_an5(i, 1) - eps)
            sign_areadiff(i, 2) = -1;
        elseif dists_f_g_an4(i, 1) > (dists_f_g_an5(i, 1) + eps)
            sign_areadiff(i, 2) = 1;
        else
            sign_areadiff(i, 2) = 0;
        end
    else
        if dists_f_g_an4(i, 1) < (dists_f_g_an5(i, 1) - eps)
            sign_areadiff(i, 2) = 1;
        elseif dists_f_g_an4(i, 1) > (dists_f_g_an5(i, 1) + eps)
            sign_areadiff(i, 2) = -1;
        else
            sign_areadiff(i, 2) = 0;
        end
    end
end


%% ############################## Calculate sign from face_itg_curv ##############################
eps_face_itg_curv = 0.1;
sign_face_itg_curv = zeros(length(tracked_uniqueface_an4), 1);
towards_left = zeros(length(tracked_uniqueface_an4), 1);
left_convex = zeros(length(tracked_uniqueface_an4), 1);
for i = 1:length(tracked_uniqueface_an4)
    mask_1 = (tracked_uniqueface_an4(i, 1) == data_face_an4(:,1) & tracked_uniqueface_an4(i, 2) == data_face_an4(:,2));
    mask_2 = (tracked_uniqueface_an4(i, 2) == data_face_an4(:,1) & tracked_uniqueface_an4(i, 1) == data_face_an4(:,2));
    if sum(mask_1) > 0
        face_itg_curv_1 = data_face_an4(mask_1, 3)*data_face_an4(mask_1, 4);
    else
        face_itg_curv_1 = 0;
    end
    if sum(mask_2) > 0
        face_itg_curv_2 = data_face_an4(mask_2, 3)*data_face_an4(mask_2, 4);
    else
        face_itg_curv_2 = 0;
    end
    
    if tracked_uniqueface_an4(i, 1) == data_face_an4(mask_1, 1)
        face_itg_curv = - face_itg_curv_1 + face_itg_curv_2;
    else
        face_itg_curv = face_itg_curv_1 - face_itg_curv_2;
    end
    
%     if face_itg_curv < - eps_face_itg_curv
%         sign_face_itg_curv(i) = 1;
%     elseif face_itg_curv > eps_face_itg_curv
%         sign_face_itg_curv(i) = -1;
%     else
%         sign_face_itg_curv(i) = 0;
%     end

    if dists_f_g_an4(i, 1) < (dists_f_g_an5(i, 1) - eps)
        towards_left(i) = -1;
    elseif dists_f_g_an4(i, 1) > (dists_f_g_an5(i, 1) + eps)
        towards_left(i) = 1;
    else
        towards_left(i) = 0;
    end
    
    if face_itg_curv > eps_face_itg_curv
        left_convex(i) = -1;
    elseif face_itg_curv < - eps_face_itg_curv
        left_convex(i) = 1;
    else
        left_convex(i) = 0;
    end
        
end



%% ############################ Write csv For Python ############################
total_dists_an4 = sum(dists_f_g_an4, 2);
total_dists_an5 = sum(dists_f_g_an5, 2);
dists_f_g_an4 = dists_f_g_an4 ./ total_dists_an4;
dists_f_g_an5 = dists_f_g_an5 ./ total_dists_an5;

% --- move_right = distance to left_grain_centroid increase ---
% move_right = (dists_f_g_an4(:,1) > (dists_f_g_an5(:,1) + eps));
% move_left = (dists_f_g_an4(:,1) < (dists_f_g_an5(:,1) - eps));
% motion = zeros(size(move_right));
% motion(move_right) = 1;
% motion(move_left) = -1;
dist_centroid_diff_l_an4 = dists_f_g_an5(:,1) - dists_f_g_an4(:,1);

% grain_size = grain_data{1,1};
% grain_nf = grain_data{1,2};
% grain_f_fnn = grain_data{1,4};

% gsize_incre_l_m_r = zeros(size(move_right));
% gf_incre_l_m_r = zeros(size(move_right));
% f_fnn_incre_l_m_r = zeros(size(move_right));
gsize_diff_an4 = zeros(size(move_right));
gf_diff_an4 = zeros(size(move_right));
f_fnn_diff_an4 = zeros(size(move_right));


for i = 1:size(tracked_uniqueface_an4, 1)
    l_idx_an4 = tracked_uniqueface_an4(i, 1);
    r_idx_an4 = tracked_uniqueface_an4(i, 2);
%     l_idx_an5 = tracked_uniqueface_an5(i, 1);
%     r_idx_an5 = tracked_uniqueface_an5(i, 2);
    
%     gsize_incre_l_m_r(i) = (grain_data_an5{1,1}(l_idx_an5) - grain_data_an4{1,1}(l_idx_an4)) - ...
%                             (grain_data_an5{1,1}(r_idx_an5) - grain_data_an4{1,1}(r_idx_an4));
%     gf_incre_l_m_r(i) = (grain_data_an5{1,2}(l_idx_an5) - grain_data_an4{1,2}(l_idx_an4)) - ...
%                             (grain_data_an5{1,2}(r_idx_an5) - grain_data_an4{1,2}(r_idx_an4));
%     f_fnn_incre_l_m_r(i) = (grain_data_an5{1,4}(l_idx_an5) - grain_data_an4{1,4}(l_idx_an4)) - ...
%                             (grain_data_an5{1,4}(r_idx_an5) - grain_data_an4{1,4}(r_idx_an4));

    gsize_diff_an4(i) = grain_data_an4{1,1}(l_idx_an4) -  grain_data_an4{1,1}(r_idx_an4);
    gf_diff_an4(i) = grain_data_an4{1,2}(l_idx_an4) -  grain_data_an4{1,2}(r_idx_an4);
    f_fnn_diff_an4(i) = grain_data_an4{1,4}(l_idx_an4) -  grain_data_an4{1,4}(r_idx_an4);

end


%%
fileID = fopen('190301_mig_signs.txt','w');
fprintf(fileID,'%s,%s,%s,%s,%s\n','gsize_diff_an4', 'gf_diff_an4', 'f_fnn_diff_an4', 'fnnf_maxadec_l_m_r', 'dist_centroid_diff_l_an4');
for i = 1:length(gsize_diff_an4)
    fprintf(fileID, '%6.3f,%6d,%6.3f,%6.3f,%6.3f\n', gsize_diff_an4(i), ...
        gf_diff_an4(i), f_fnn_diff_an4(i), fnnf_maxadec_l_m_r(i), dist_centroid_diff_l_an4(i));
end
fclose(fileID);





