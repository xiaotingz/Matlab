function plotTriNormals(normal_svm, obj_facelabel_an4, obj_facelabel_an5, file_an4, file_an5, obj)
% ############################################################################
% Input
%   - obj_face_an, [2, 1] the objective face label
%       doesn't need to be sorted
%   - obj, string, the choice of showing the average plane or not
%     'show_avg'
% Notes
%   - use together with visualizeFace.m and plotSVMPlane.m
% ############################################################################
% ------------------ load data for debug --------------------
% obj_face_an4 = face_node_info{1,1};
% obj_face_an5 = face_node_info{1,2};
% -----------------------------------------------------------
color1 = [0, 0.4470, 0.7410];
color2 = [0.9290, 0.6940, 0.1250];
color3 = [0.4660, 0.6740, 0.1880];


    % ----- Read In Data ----- 
    normal_an4 = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals')';
    normal_an5 = h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals')';
    centroid_an4 = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceCentroids')';
    centroid_an5 = h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceCentroids')';
    facelabel_an4 = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
    facelabel_an5 = h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';

    mask_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2)) ...
     | (facelabel_an4(:,1) == obj_facelabel_an4(2) & facelabel_an4(:,2) == obj_facelabel_an4(1));
    mask_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2)) ... 
     | (facelabel_an5(:,1) == obj_facelabel_an5(2) & facelabel_an5(:,2) == obj_facelabel_an5(1));

    centroid_an4_plot = centroid_an4(mask_an4, :);
    centroid_an5_plot = centroid_an5(mask_an5, :);    
    normal_an4_plot = normal_an4(mask_an4, :);
    normal_an5_plot = normal_an5(mask_an5, :);
    
    % ----- Change winding and Recombine to plot -----
    mask_reverse_an4 = facelabel_an4(mask_an4,1) > facelabel_an4(mask_an4,2);
    mask_reverse_an5 = facelabel_an5(mask_an5,1) > facelabel_an5(mask_an5,2);
    normal_an4_plot(mask_reverse_an4, :) = - normal_an4_plot(mask_reverse_an4, :);
    normal_an5_plot(mask_reverse_an5, :) = - normal_an5_plot(mask_reverse_an5, :);

    % ----- Plot triangle normals-----
    quiver3(centroid_an4_plot(:,1), centroid_an4_plot(:,2), centroid_an4_plot(:,3), normal_an4_plot(:,1), ...
        normal_an4_plot(:,2), normal_an4_plot(:,3), 'color', color1,  'LineWidth', 1, 'MaxHeadSize', 3);
    hold on
    rotate3d on
    quiver3(centroid_an5_plot(:,1), centroid_an5_plot(:,2), centroid_an5_plot(:,3), normal_an5_plot(:,1), ...
        normal_an5_plot(:,2), normal_an5_plot(:,3), 'color', color2,  'LineWidth', 1, 'MaxHeadSize', 3);
    
    % ----- Prepare for the grid -----
    all_centroids = [centroid_an4_plot; centroid_an5_plot];
    point = sum(all_centroids) / length(all_centroids);
    x_min = min(all_centroids(:,1));
    x_max = max(all_centroids(:,1));
    y_min = min(all_centroids(:,2));
    y_max = max(all_centroids(:,2));
    [xx,yy] = ndgrid(x_min : (x_max - x_min)/2 : x_max, y_min : (y_max - y_min)/2 : y_max);
    
    normal_an4_avg = sum(normal_an4_plot)/length(normal_an4_plot);
    normal_an5_avg = sum(normal_an5_plot)/length(normal_an5_plot);
    % ----- a plane is a*x+b*y+c*z+d=0, [a,b,c] is the normal. ----- 
    % ----- Thus, we have to calculate d and we're set ----- 
    d_an4 = - dot(point, normal_an4_avg);
    d_an5 = - dot(point, normal_an5_avg);
    
    z_an4 = (-normal_an4_avg(1) * xx - normal_an4_avg(2) * yy - d_an4) / normal_an4_avg(3);
    z_an5 = (-normal_an5_avg(1) * xx - normal_an5_avg(2) * yy - d_an5) / normal_an5_avg(3);
    
    d_svm = - dot(point, normal_svm);
    z_svm = (-normal_svm(1) * xx - normal_svm(2) * yy - d_svm) / normal_svm(3);
    scale = ceil(max(length(normal_an4_plot), length(normal_an5_plot))/50);
    scale(scale < 8) = 8;
    
    idx = 2;
    quiver3(xx(idx)+0.1, yy(idx), z_svm(idx), normal_svm(1), normal_svm(2), normal_svm(3), scale, 'color', color3,  'LineWidth', 3, 'MaxHeadSize', 3);
    quiver3(xx(idx)+0.1, yy(idx), z_svm(idx), -normal_svm(1), -normal_svm(2), -normal_svm(3), scale, 'color', color3, 'LineWidth', 3,'MaxHeadSize', 3);
    surf(xx,yy,z_svm, 'FaceColor', color3, 'EdgeColor', color3, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);

if nargin == 6 && strcmp(obj, 'show_avg')
    % ----- plot the plane corresponding to the average plane normal -----
    surf(xx,yy,z_an4, 'FaceColor', color1, 'EdgeColor', color1, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
    surf(xx,yy,z_an5, 'FaceColor', color2, 'EdgeColor', color2, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
    % ----- plot the average plane normal -----
%     scale = ceil(max(length(normal_an4_plot), length(normal_an5_plot))/50);
%     scale(scale < 8) = 8;
%     idx = 5;
%     quiver3(xx(idx), yy(idx), z_an4(idx), normal_an4_avg(1), normal_an4_avg(2), normal_an4_avg(3), scale, 'color', color1,  'LineWidth', 3, 'MaxHeadSize', 3);
%     quiver3(xx(idx), yy(idx), z_an5(idx), normal_an5_avg(1), normal_an5_avg(2), normal_an5_avg(3), scale, 'color', color2, 'LineWidth', 3,'MaxHeadSize', 3);
    % ----- set plot ax aspect ratio -----
    daspect([1 1 1])
%     plotSVMPlane(normal_an4_avg, d_an4, [x_min, x_max], [y_min, y_max], scale)
end    

end
