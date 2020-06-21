function [face_dist_avgs, face_dist_stds, tri_fls, tri_dists] = calcDistFromTargetNormals(file, obj_faces, target_normals, eps_curv, eps_area, eps_min_ang)
% ###########################################################################
% * Input
%     - obj_faces = [n,2] = [label_a, label_b]
%         Facelabels of the faces of interest
%     - target_normals = [m, 3]
%         Plane normals of interest, e.g. [1,1,1], [1,1,0], [1,0,0]
%     - eps_
%         Filter bad triangles, this requires corresponding D3D filters to be run.
% * Output
%     - face_dist_avgs = [n,1] 
%         Average distance of triangle normals from target_normals
%     - face_dist_stds = [n, 1]
%         Standard deviation from the distances.
%     - tri_fls = [k, 2]
%         Face normal of the triangles
%     - tri_dists = [k, 3]
%         Distance from the triangles to the target normals.
%         Output just for debug purpose.
% * Note 
%     - This file is to be used in featurePrepOther.m.
%     - Dependency: calcTriDistFromTargetNormals.m
%     - This function aims to characterize the plane normal of a grain
%     face in a few scalar numbers. Scalars because we don't want deep NN.
%     A grain face include many different normal directions, so the scalars
%     are chosen as the avg(dist(triangle_normals, target_normals)) and 
%     std(dist(triangle_normals, target_normals)).
% ###########################################################################
% % ------------------ load data for debug --------------------
% file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
% eps_curv = 1;
% eps_area = 7;
% eps_min_ang = 10;
% % load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/trials/data_matlab/190624_tracked_faces_full.mat', 'tracked_uniqueface_an4_full')
% % faces = tracked_uniqueface_an4_full;
% % clear tracked_uniqueface_an4_full
% target_normals = [1,1,1;1,1,0;1,0,0;1,1,2];
% target_normals = target_normals ./ vecnorm(target_normals, 2, 2);
% % -----------------------------------------------------------

% ----- Read data -----
tri_normals = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/FaceNormals').';
tri_fls = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels').';
ea =  h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles').';
ea(1, :) = [];
ea = rad2deg(ea);

% ----- Threshold bad and triple line triangles -----
tri_nodes = 1 + h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')';
node_types = h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType')';
tri_curv = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
tri_area = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
tri_min_ang = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
tri_node_types = node_types(tri_nodes);
mask_good = all(tri_fls > 0, 2) & all(tri_node_types == 2, 2) & ... 
    abs(tri_curv) < eps_curv & tri_area < eps_area & tri_min_ang > eps_min_ang;
clear tri_area tri_curv tri_min_ang tri_node_types tri_nodes 

tri_fls = tri_fls(mask_good, :);
tri_normals = tri_normals(mask_good, :);

% ----- For every triangle, calculate distance to the target normals -----
tri_dists = calcTriDistFromTargetNormals_formal(tri_normals, tri_fls, ea, target_normals);
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Lsmoth_normal_dists.mat')

% ----- For every face, calculate average distance and standard devation -----

tri_fls = sort(tri_fls, 2);
obj_faces = sort(obj_faces, 2);
face_dist_avgs = ones(size(obj_faces, 1), size(target_normals, 1)) * NaN;
face_dist_stds = ones(size(obj_faces, 1), size(target_normals, 1)) * NaN;
for i = 1:size(obj_faces, 1)
    mask = (tri_fls(:, 1)==obj_faces(i, 1) & tri_fls(:, 2)==obj_faces(i, 2));
    tri_on_face = tri_dists(mask, :);
    
    face_dist_avgs(i, :) = sum(tri_on_face, 1) / size(tri_on_face, 1);
    % """ std(data, weight_factor, dim), default weight_factor = 0 """    
    face_dist_stds(i, :) = std(tri_on_face, 0, 1);

end

end

% % ###################### Plots ######################
% h1 = cdfplot(face_dist_avgs(:,1));
% hold on
% h2 = cdfplot(face_dist_avgs(:,2));
% h3 = cdfplot(face_dist_avgs(:,3));
% h4 = cdfplot(face_dist_avgs(:,4));
% title('avg\_dist from target normal')
% xlabel('average distance, °')
% 
% set(h1, 'LineWidth', 2);
% set(h2, 'LineWidth', 2);
% set(h3, 'LineWidth', 2);
% set(h4, 'LineWidth', 2);
% 
% legend({'[111]', '[110]', '[100]', '[112]'}, 'Location', 'SouthEast')
% print('dist_from_target_normals', '-dtiff','-r300')
% 
% scatter(face_dist_avgs(:,1), face_dist_stds(:,1), 'filled', 'MarkerFaceAlpha', 0.3)
% ylabel('std, °')
% xlabel('average distance from [111], °')
% print('std_dist_[111]', '-dtiff','-r300')
% 
% scatter(face_dist_avgs(:,2), face_dist_stds(:,2), 'filled', 'MarkerFaceAlpha', 0.3)
% ylabel('std, °')
% xlabel('average distance from [110], °')
% print('std_dist_[110]', '-dtiff','-r300')




