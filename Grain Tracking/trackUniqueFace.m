% function [tracked_uniqueface_1, tracked_uniqueface_2] = trackUniqueFace(file_1, file_2, look_up_table, use_complete_faces)
% % ###########################################################################
% % * Very similar to TrackFace, except that the input is the uniqueFaces.
% % * Input
% %     - lookUp = [N,2] = [ID_s1, ID_s2] 
% %         is a table which contains only the tracked faces and was sorted in the order of state2.
% %     - complete, logical variable
% %         =0 --> all faces; 
% %         =1 --> the tracked faces are complete in both anneal states
% % * Output
% %     - tracked_uniqueface_1 and tracked_uniqueface_2 = [m, 2]
% %         m is the number of unique facesd. The two arrays have been sorted
% %         such that every row corresponds to the same face. 
% %         this correspondence can be applied to faceCoords directly.
% % ###########################################################################
% % ------------------ load data for debug --------------------
% % file_1 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
% % file_2 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
% % load('look_up_table_an4_an5.mat');
% % -----------------------------------------------------------
%     unique_facelabel_1 = h5read(file_1, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
%     unique_facelabel_2 = h5read(file_2, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
%     unique_facelabel_1 = unique_facelabel_1(all(unique_facelabel_1>0, 2), :);
%     unique_facelabel_2 = unique_facelabel_2(all(unique_facelabel_2>0, 2), :);
%     unique_facelabel_1 = sortrows(unique_facelabel_1);
%     unique_facelabel_2 = sortrows(unique_facelabel_2);
%     if strcmp(use_complete_faces, 'use_complete_faces')
%             surface_grain_1 = h5read(file_1,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
%             surface_grain_2 = h5read(file_2,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
%             surface_grain_1(1) = [];
%             surface_grain_2(1) = [];
%             % --- The face is incomplete if both grains are on free surface ---
%             mask_completeface_1 = (sum(surface_grain_1(unique_facelabel_1), 2) < 2);
%             mask_completeface_2 = (sum(surface_grain_2(unique_facelabel_2), 2) < 2);
%             unique_facelabel_1 = unique_facelabel_1(mask_completeface_1, :);
%             unique_facelabel_2 = unique_facelabel_2(mask_completeface_2, :);    
%     else
%         warning('All faces are tracked. Add in a fourth argument if want Complete Faces only.');
%     end
% 
% % ##### Make a Full Look-up Table ##### 
%     look_up_table = sortrows(look_up_table, 2);  %  somehow the look up table is not really sorted well 
%     look_up_2to1 = zeros(max(look_up_table(:,2)),1);
%     idx = 1;
%     for i = 1 : max(look_up_table(:,2))
%         if look_up_table(idx,2) == i
%             look_up_2to1(i,1) = look_up_table(idx,1);
%             idx = idx + 1;
%         else
%             look_up_2to1(i,1) = NaN;
%         end
%     end
%     
%     unique_facelabel_1 = [(1:length(unique_facelabel_1)).', unique_facelabel_1];
%     unique_facelabel_2in1 = [(1:length(unique_facelabel_2)).', look_up_2to1(unique_facelabel_2)];
%     mask = ~(any(isnan(unique_facelabel_2in1),2));
%     unique_facelabel_2in1 = unique_facelabel_2in1(mask,:);
%     unique_facelabel_2in1(:,2:3) = sort(unique_facelabel_2in1(:,2:3),2);
%     unique_facelabel_2in1 = sortrows(unique_facelabel_2in1,[2,3]);
% 
% % ##### Track Faces by Comparing Grain IDs #####
%     length_1 = length(unique_facelabel_1);
%     length_2 = length(unique_facelabel_2in1);
%     i = 1;
%     j = 1;
%     unique_facelabel_corresp = [];
%     while j <= length_1 && i <= length_2
%         if unique_facelabel_2in1(i,2) == unique_facelabel_1(j,2)
%             if unique_facelabel_2in1(i,3) == unique_facelabel_1(j,3)
%                 unique_facelabel_corresp = vertcat(unique_facelabel_corresp, [unique_facelabel_1(j,1), unique_facelabel_2in1(i,1)]);
%                 i = i + 1;
%                 j = j + 1;
%             elseif unique_facelabel_2in1(i,3) > unique_facelabel_1(j,3)
%                 j = j + 1;
%             else
%                 i = i + 1;
%             end
%         elseif unique_facelabel_2in1(i,2) > unique_facelabel_1(j,2)
%             j = j + 1;
%         else
%             i = i + 1;
%         end
%     end
%     unique_facelabel_corresp = sortrows(unique_facelabel_corresp, 2);
%     
%     tracked_uniqueface_1 = unique_facelabel_1(unique_facelabel_corresp(:,1), 2:3);
%     tracked_uniqueface_2 = unique_facelabel_2(unique_facelabel_corresp(:,2), :);
% end
% 
% 
% 

% % ###################################### Check, Or Method2 ######################################
file_1 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin3_Hsmooth.dream3d';
file_2 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4_Hsmooth.dream3d';
load('look_up_table_an4_an5crop.mat')

% ------------- read in facelabels -------------
fl_1 = h5read(file_1, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
fl_2 = h5read(file_2, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
surface_grain_1 = logical(h5read(file_1,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')).';
surface_grain_2 = logical(h5read(file_2,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')).';
surface_grain_1(1) = [];
surface_grain_2(1) = [];
ids_1 = (1:length(surface_grain_1))';
ids_2 = (1:length(surface_grain_2))';
surface_grain_1 = ids_1(surface_grain_1);
surface_grain_2 = ids_2(surface_grain_2);

% ------------- filter out bad and surface faces -------------
mask_surf_face_1 = all(ismember(fl_1, surface_grain_1), 2);
mask_surf_face_2 = all(ismember(fl_2, surface_grain_2), 2);
mask_face_1 = all(fl_1 > 0, 2) & ~ mask_surf_face_1;
mask_face_2 = all(fl_2 > 0, 2) & ~ mask_surf_face_2;
% mask_face_1 = all(fl_1 > 0, 2);
% mask_face_2 = all(fl_2 > 0, 2);

fl_1 = fl_1(mask_face_1, :);
fl_2 = fl_2(mask_face_2, :);
fl_1 = sort(fl_1, 2);
fl_2 = sort(fl_2, 2);

% ------------- construct fl_2from1 from fl_11 -------------
map_1to2 = containers.Map(look_up_table(:,1), look_up_table(:,2));
fl_2from1 = - ones(size(fl_1));
for i = 1:size(fl_1, 1)
    for j = 1:2
        if isKey(map_1to2, fl_1(i, j))
            fl_2from1(i, j) = map_1to2(fl_1(i, j));
        end
    end
end

% ------------- sort fl_2from1 and determine a mask by comparing fl_2from1 and fl_2 -------------
fl_2from1 = sort(fl_2from1, 2);
mask = ismember(fl_2from1, fl_2, 'rows');

% ------------- return fl_1(mask, :) and fl_2from1(mask, :) -------------
tracked_uniqueface_1 = fl_1(mask, :);
tracked_uniqueface_2 = fl_2from1(mask, :);

% sum(ismember(tracked_uniqueface_1, tracked_uniqueface_an4, 'rows'))
% sum(ismember(tracked_uniqueface_2, tracked_uniqueface_an5, 'rows'))



