file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4.dream3d';
load('../Grain Tracking/look_up_table_an4_an5crop.mat')
% file_an4 = '/Volumes/XIAOTING/Ni/simu_An4_clean_seg.dream3d';
% file_an5 = '/Volumes/XIAOTING/Ni/simu_An5_clean_seg.dream3d';

eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;

faces_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
faces_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
faces_an4(1, :) = [];
faces_an5(1, :) = [];
faces_an4 = faces_an4(all(faces_an4>0, 2), :);
faces_an5 = faces_an5(all(faces_an5>0, 2), :);
faces_an4 = sortrows(sort(faces_an4, 2));
faces_an5 = sortrows(sort(faces_an5, 2));


[result_an4, ~, ~, ~] = findQuadNodes(file_an4);
[triple_line_an4, tl_info_an4] = findTripleLines(file_an4, eps_area, eps_curv, eps_min_ang);
[triple_line_full_an4, tl_info_full_an4] = findTripleLines(file_an4, 10e6, 10e6, 0);
[num_corners_an4, num_edges_an4] = getFaceCharacter(faces_an4, triple_line_full_an4, result_an4{1}, result_an4{2});
% [num_nnface_avgcorner_an4] = findFaceNNAvgCorner(file_an4, faces_an4, num_corners_an4, triple_line_full_an4);

[result_an5, ~, ~, ~] = findQuadNodes(file_an5);
[triple_line_an5, tl_info_an5] = findTripleLines(file_an5, eps_area, eps_curv, eps_min_ang);
[triple_line_full_an5, tl_info_full_an5] = findTripleLines(file_an5, 10e6, 10e6, 0);
[num_corners_an5, num_edges_an5] = getFaceCharacter(faces_an5, triple_line_full_an5, result_an5{1}, result_an5{2});
% [num_nnface_avgcorner_an5] = findFaceNNAvgCorner(file_an5, faces_an5, num_corners_an5, triple_line_full_an5);

% [avg_nng_diff, max_nng_diff] = findFaceLocalTopologyChange(file_an4, file_an5, faces_an4, look_up_table, triple_line_full_an4);


% [da_len_weighted_an4, da_num_weighted_an4] = calcGrainFaceDAs(faces_an4, triple_line_an4, tl_info_an4);
% [da_len_weighted_an5, da_num_weighted_an5] = calcGrainFaceDAs(faces_an5, triple_line_an5, tl_info_an5);

%% ##################################### Checks #####################################
% Most of the checks are at the end of each function file.
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/final/190730_Hsmoth_GeoTopo.mat', ...
    'face_area_an4', 'face_area_diff', 'face_corners_an4', 'face_edges_an4', 'tracked_uniqueface_an4_inner');
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/final/190730_Hsmoth_MeanField.mat', ...
    'tracked_uniqueface_an4_complete')
% ('face_itg_abscurv_an4', 'face_itgcurv_an4_left')

face_area_an5 = face_area_an4 + face_area_diff;
area_diff_ratio = face_area_diff ./ face_area_an4;
mask_good = face_area_an4 > 20 & face_area_an5 > 20 & area_diff_ratio < 10 & area_diff_ratio > -0.9;
mask_complete = ismember(tracked_uniqueface_an4_inner, tracked_uniqueface_an4_complete, 'rows');
mask = mask_good & mask_complete;
face_area_an4 = face_area_an4(mask);
face_edges_an4 = face_edges_an4(mask);
face_corners_an4 = face_corners_an4(mask);

clear mask mask_complete mask_good tracked_uniqueface_an4_complete tracked_uniqueface_an4_inner face_area_an5 face_area_diff area_diff_ratio

%% ###################### #edges v.s. #corners ######################
set(0,'defaultAxesFontSize',20)
% scatter(face_corners_an4, face_edges_an4, 'filled', ...
%     'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.1)
% xlabel('C')
% ylabel('E')
% axis_lim = floor(min(max(face_corners_an4), max(face_edges_an4)) / 5 + 1) * 5;
% line([0, axis_lim], [0, axis_lim], 'color', [0.5, 0.5, 0.5])
% xlim([0, axis_lim]);
% ylim([0, axis_lim]);
% pbaspect([1 1 1])

% figure()
% x = face_corners_an4;
% y = face_edges_an4;
% x_bin_setting = 0 : axis_lim : axis_lim;
% y_bin_setting = 0 : axis_lim : axis_lim;
% plotDensityMatrix(x, y, x_bin_setting, y_bin_setting, 'C', 'E')
% 
% xTicks = [0:5:35];
% % xLabels = mat2cell();
% set(gca, 'XTick', xTicks)
% set(gca,'XMinorTick','on','YMinorTick','on')


figure()
x = sqrt(face_area_an4);
y = face_edges_an4;
step_x = 6;
edges_x = 0 : step_x : step_x*axis_lim;
edges_y = 0 : 1 : axis_lim;
plotDensityMatrix(x, y, edges_x, edges_y, 'sqrt(A)', 'E')


xlabel('$\sqrt{A}$  ($\mu$m)', 'Interpreter', 'Latex', 'FontSize',20);
% xTicks = [0:5:35];
% xLabels = mat2cell();
% set(gca, 'XTick', xTicks)
xTickLabels = [5:5:35]*step_x;
xTickLabels = num2cell(xTickLabels,length(xTickLabels));
xticklabels(xTickLabels)

set(gca,'XMinorTick','on','YMinorTick','on')



%%





