% function grain_area = calcGrainArea(data_face)
% ##########################################################################
% * Input
%     - data_face = [m,4]  
%           [g1, g2, area, avg_cur], returned by calling calcFaceCurvature.m in G_F_mF.m in Grain Curvature
% * Output
%     - areas = [n, 2]
%           grain areas, corresponding to increasing id order
% * Note
%     - related functions: calcFaceToGrainCentroidDist.m
% ##########################################################################
% ----------------------- load debug data -----------------------
% load('190228_nn_faces.mat');
data_face = data_face_an5;
% ---------------------------------------------------------------

% ##### Calcualte Area for All Related Grains ##### 
grains = unique(data_face(:,1:2));
areas = zeros(size(grains));
for i = 1:length(grains)
    mask = (data_face(:,1) == grains(i) | data_face(:,2) == grains(i));
    face_areas = data_face(mask, 3);
    areas(i) = sum(face_areas);
end
% 
% area_map = containers.Map(grains, areas);
% 
% 
% % ##### Order Data for Output ##### 
% grain_area = zeros(size(grains));
% for i = 1:size(grains,1)
%     
%         grain_area(i, j) = area_map(faces(i, j));
%     end
% end









