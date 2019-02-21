% function mig_signs = calcFaceMigSign(file_an4, file_an5, tracked_uniqueface_an4, tracked_uniqueface_an5, energy_an4, energy_an5)
% ##########################################################################
% * Input
%     - tracked_uniqueface_ = [n,2]  
%           returned by trackFace.m or trackUniqueFace.m
%     - energy_ = {[m, 1], k}
%           m = #grains, k = #energies passed in
%           size, or F - <Fnn>, or integral curvature.
% * Output
%     - mig_signs = [n, k]
%           n = #grain_faces, k = #energies passed in
% * NOTE
%     - This function is to calculate the sign of face migration. 
%     - sign = 1    if face moving towards the grain with lower energy; 
%       sign = -1   if face moving toward the grain with higher energy;
%       sign = 0    if the two grains swapped energy level.
%     - Dependence: calcFaceToGrainCentroidDist.m
% ##########################################################################
% ----------------------- load debug data -----------------------
% file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
% file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/181107.mat', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5');
% load('look_up_table_an4_an5.mat')

% ##### run  /Grain Curvature/G_F_mF.m to get the following data: from data_grain & F_mF_diff #####
% energy_an4 = {grain_size_an4, grain_F_an4, grain_itg_curv_an4, F_avgFnn_an4};
% energy_an5 = {grain_size_an5, grain_F_an5, grain_itg_curv_an5, F_avgFnn_an5};
% ---------------------------------------------------------------
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

% ##### Calc Face To Grain Centroid Dist #####
% ----- Absolute distance -----
% dists_f_g_an4 = calcFaceToGrainCentroidDist(file_an4, tracked_uniqueface_an4);
% dists_f_g_an5 = calcFaceToGrainCentroidDist(file_an5, tracked_uniqueface_an5);
% ----- Convert distance to portion -----
total_dists_an4 = sum(dists_f_g_an4, 2);
total_dists_an5 = sum(dists_f_g_an5, 2);
dists_f_g_an4 = dists_f_g_an4 ./ total_dists_an4;
dists_f_g_an5 = dists_f_g_an5 ./ total_dists_an5;

% % ##### Check Sign For Each Grain Energy Metric #####
eps = 1e-2;
sign = zeros(size(tracked_uniqueface_an4, 1), size(energy_an4, 2));
i = 1;
for j = 1:size(tracked_uniqueface_an4, 1)
    energy_an4_cur = energy_an4{i};
    energy_an5_cur = energy_an5{i};
    id_l_an4 = tracked_uniqueface_an4(j, 1);
    id_r_an4 = tracked_uniqueface_an4(j, 2);
    id_l_an5 = tracked_uniqueface_an5(j, 1);
    id_r_an5 = tracked_uniqueface_an5(j, 2);
    % """
    % eps: a small variable to filter out numerical errors.
    % smaller size = bigger energy ??? 
    % only energy in an4 should be considered, or energy in an4 & an5 both
    % needs to be considered?
    % """
%     left_grain_energy_big = (energy_an4_cur(id_l_an4) < energy_an4_cur(id_r_an4) && ...
%                             energy_an5_cur(id_l_an5) < energy_an5_cur(id_r_an5));
%     right_grain_energy_big = (energy_an4_cur(id_l_an4) > energy_an4_cur(id_r_an4) && ...
%                             energy_an5_cur(id_l_an5) > energy_an5_cur(id_r_an5));
    left_grain_small = (energy_an4_cur(id_l_an4) < energy_an4_cur(id_r_an4));
    if left_grain_small
        if dists_f_g_an4(j, 1) < (dists_f_g_an5(j, 1) - eps)
            sign(j, i) = 1;
        elseif dists_f_g_an4(j, 1) > (dists_f_g_an5(j, 1) + eps)
            sign(j, i) = -1;
        else
            sign(j, i) = 0;
        end
    else
        if dists_f_g_an4(j, 1) < (dists_f_g_an5(j, 1) - eps)
            sign(j, i) = -1;
        elseif dists_f_g_an4(j, 1) > (dists_f_g_an5(j, 1) + eps)
            sign(j, i) = 1;
        else
            sign(j, i) = 0;
        end
    end

end









