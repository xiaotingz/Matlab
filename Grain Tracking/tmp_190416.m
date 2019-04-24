%% ###################### Unify Curvature Direction ##########################
file = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin2_smooth_curveswaped.dream3d';
curv = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures')';
fl = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
mask_inner = all(fl > 0, 2);
fl_inner = fl(mask_inner, :);
curv_inner = curv(mask_inner, :);

% """ 
% if obj_grain on left of face_label, tri_curve = - tri_curve
% """
for i = 1:size(tracked_uniqueface_an4, 1)
    disp(i)
%     mask = ismember(fl, tracked_uniqueface_an4(i,:), 'rows');
    mask = (fl_inner(:,1) == tracked_uniqueface_an4(i,1) & fl_inner(:,2) == tracked_uniqueface_an4(i,2));
    curv_inner(mask) = - curv_inner(mask);
end
    
curv(mask_inner) = curv_inner;
h5write(file, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures', curv');

%% ###################### Calculate face centroid of interest ##########################
load('/Volumes/XIAOTING/Ni/working/Grain Tracking/tmp_190416.mat', ...
    'file_an4', 'file_an5', 'tracked_uniqueface_an4_sort', 'tracked_uniqueface_an5');
% ##### Read and Clean Data #####
fl_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
tri_node_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
fl_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
tri_node_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

tri_node_an4 = tri_node_an4(all(fl_an4 > 0, 2), :);
fl_an4 = fl_an4(all(fl_an4 > 0, 2), :);
fl_an4 = sort(fl_an4, 2);
tri_node_an5 = tri_node_an5(all(fl_an5 > 0, 2), :);
fl_an5 = fl_an5(all(fl_an5 > 0, 2), :);
fl_an5 = sort(fl_an5, 2);

% ##### Id for faces of interest #####
obj_idx = [4484, 5623];
obj_face_an4 = tracked_uniqueface_an4_sort(obj_idx, :);
obj_face_an5 = tracked_uniqueface_an5(obj_idx, :);

% ##### Find face_centroid #####
face_centroid_an4 = zeros(size(obj_face_an4, 1), 3);
face_centroid_an5 = zeros(size(obj_face_an4, 1), 3);
for i = 1:size(obj_face_an4)
    mask_an4 = (fl_an4(:,1) == obj_face_an4(i, 1) & fl_an4(:,2) == obj_face_an4(i, 2) ...
             | fl_an4(:,1) == obj_face_an4(i, 2) & fl_an4(:,2) == obj_face_an4(i, 1));
    nodes_an4 = unique(tri_node_an4(mask_an4, :));
    face_centroid_an4(i, :) = sum(node_coord_an4(nodes_an4, :))/length(nodes_an4);
    mask_an5 = (fl_an5(:,1) == obj_face_an5(i, 1) & fl_an5(:,2) == obj_face_an5(i, 2) ...
             | fl_an5(:,1) == obj_face_an5(i, 2) & fl_an5(:,2) == obj_face_an5(i, 1));
    nodes_an5 = unique(tri_node_an5(mask_an5, :));
    face_centroid_an5(i, :) = sum(node_coord_an5(nodes_an5, :))/length(nodes_an5);
end

diff = face_centroid_an4 - face_centroid_an5;





load('/Volumes/XIAOTING/Ni/working/Grain Tracking/tmp_190416.mat', ...
    'file_an4', 'file_an5', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5');

rfvec_twin = ([1,1,1]/norm([1,1,1]))*tand(60/2);
eps = 0.05;
mask = vecnorm(rfvecs_an4 - rfvec_twin, 2, 2)<eps & vecnorm(rfvecs_an5 - rfvec_twin, 2, 2)<eps;

idx_grains = (1:size(dist_ctwin, 1))';
idx_twins = idx_grains(mask);






%% ###################### Plot Twins ##########################
rng('shuffle');
i = randi(length(idx_twins));
idx = idx_twins(i);


face_node_info = getSingleFaceNodes(tracked_uniqueface_an4(idx,:), tracked_uniqueface_an5(idx,:));

figure
visualizeFace(face_node_info)
disp(['pair ', num2str(idx), ',   dist_ctwin=', num2str(dist_ctwin(idx))]);


f_name = ['twins_pair_', num2str(idx)];
% f_name = ['twins_pair_', num2str(idx_twins(idx)), '_2'];
print(f_name, '-dtiff', '-r300')







%% ###################### Display data & Make plot ##########################
% load('/Volumes/XIAOTING/Ni/working/Grain Tracking/projection/190404_proj_dists.mat', ...
%     'face_onepiece', 'mig_localnorm_nearest', 'mig_localnorm_ot');
load('/Volumes/XIAOTING/Ni/working/Grain Tracking/projection/190417_dist_left_from_an4normproj.mat');
load('/Volumes/XIAOTING/Ni/working/Grain Tracking/tmp_190416.mat')

% daspect([1 1 1]);
% ##### plot the centroids #####
rng('shuffle');
% i = randi(size(face_onepiece, 1), 1);  
i = 1638;
idx = face_onepiece(i);

eps = 0.5;
if dist_left(idx) > eps
    memo = '   (move left)';
elseif dist_left(idx) < -eps
    memo = '   (move right)';
else
    memo = '   (stable)';
end
if gV_diff_an4(idx) > 0
    memo2 = 'right big,   ';
else
    memo2 = 'left big,   ';
end

disp(' ');
disp('#############################')
disp(['Pair ', num2str(idx), '  ;  one_piece_id = ', num2str(i)]);
dispFacePairInfo(file_an4, file_an5, tracked_uniqueface_an4, tracked_uniqueface_an5, idx)
% disp(['mig_sign_norm_nearest = ',num2str(mig_localnorm_nearest(i, 3:4))])
% disp(['mig_sign_norm_ot = ',num2str(mig_localnorm_ot(i, 3:4))])
disp('---------------')
disp(['dist_left_norm_proj = ',num2str(dist_left(idx)), memo])
disp('---------------')
disp(['mig_left_centrod = ', num2str(move_left(idx)), ...
    ',  left_centrod_an4 = ',num2str(dists_f_g_an4(idx, 1)), ',  left_centrod_an5 = ',num2str(dists_f_g_an5(idx, 1))])
disp('---------------')
disp([memo2, 'gV_diff_an4 = ', num2str(gV_diff_an4(idx)), ',  fMs_an4_left = ',num2str(fMs_an4_left(idx, 1))])
disp('---------------')
disp(['gF_diff_an4 = ', num2str(gF_diff_an4(idx)), ',  gMs_diff_an4 = ',num2str(gMs_diff_an4(idx, 1))])


figure
plotCentroids(file_an4, file_an5, tracked_uniqueface_an4, tracked_uniqueface_an5, idx)
% if move_left(idx) > 0
%     disp('move left')
%     title('move left', 'fontSize', 18)
% elseif move_left(idx) < 0
%     disp('move right')
%     title('move right', 'fontSize', 18)
% else
%     disp('stable')
%     title('stable', 'fontSize', 18)
% end

%% ##### save plot #####
f_name = ['centroids_pair_', num2str(idx)];
print(f_name, '-dtiff', '-r300')








