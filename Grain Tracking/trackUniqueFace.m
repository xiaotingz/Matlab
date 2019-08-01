function [tracked_uniqueface_1, tracked_uniqueface_2] = trackUniqueFace(file_1, file_2, look_up_table, faces_to_use)
% ###########################################################################
% * Very similar to TrackFace, except that the input is the uniqueFaces.
% * Input
%     - lookUp = [N,2] = [ID_s1, ID_s2] 
%         is a table which contains only the tracked faces and was sorted in the order of state2.
%     - faces_to_use
%         use_inner_faces --> all inner faces; 
%         use_complete_faces --> the tracked faces are complete in both anneal states
%         use_all_true_faces --> use all true faces, including the surface ones
% * Output
%     - tracked_uniqueface_1 and tracked_uniqueface_2 = [m, 2]
%         m is the number of unique facesd. The two arrays have been sorted
%         such that every row corresponds to the same face. 
%         this correspondence can be applied to faceCoords directly.
% ###########################################################################
% ------------------ load data for debug --------------------
% file_1 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
% file_2 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4.dream3d';
% load('look_up_table_an4_an5crop.mat');
% use_complete_faces = 'use_complete_faces';
% -----------------------------------------------------------
% 
% ###################################### Method1 ######################################
%     unique_facelabel_1 = h5read(file_1, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
%     unique_facelabel_2 = h5read(file_2, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
%     unique_facelabel_1 = unique_facelabel_1(all(unique_facelabel_1>0, 2), :);
%     unique_facelabel_2 = unique_facelabel_2(all(unique_facelabel_2>0, 2), :);
%     unique_facelabel_1 = sortrows(unique_facelabel_1);
%     unique_facelabel_2 = sortrows(unique_facelabel_2);
%     if strcmp(faces_to_use, 'use_complete_faces')
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
%     look_up_2to1 = zeros(max(unique_facelabel_2(:)),1);
%     idx = 1;
%     for i = 1 : max(unique_facelabel_2(:))
%         if look_up_table(idx,2) == i && idx < size(look_up_table, 1)
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



% % ###################################### Check, Or Method2 ######################################
% ------------- read in facelabels -------------
faces_1 = h5read(file_1, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
faces_2 = h5read(file_2, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';

% ------------- filter out bad and surface faces -------------
if strcmp(faces_to_use, 'use_complete_faces')
    surface_grain_1 = logical(h5read(file_1,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')).';
    surface_grain_2 = logical(h5read(file_2,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')).';
    surface_grain_1(1) = [];
    surface_grain_2(1) = [];
    ids_1 = (1:length(surface_grain_1))';
    ids_2 = (1:length(surface_grain_2))';
    surface_grain_1 = ids_1(surface_grain_1);
    surface_grain_2 = ids_2(surface_grain_2);

    mask_surf_face_1 = all(ismember(faces_1, surface_grain_1), 2);
    mask_surf_face_2 = all(ismember(faces_2, surface_grain_2), 2);
    mask_1 = all(faces_1 > 0, 2) & ~ mask_surf_face_1;
    mask_2 = all(faces_2 > 0, 2) & ~ mask_surf_face_2;
elseif strcmp(faces_to_use, 'use_all_true_faces')
    mask_1 = all(faces_1 >= 0, 2) & any(faces_1 > 0, 2);
    mask_2 = all(faces_2 >= 0, 2) & any(faces_2 > 0, 2);
    look_up_table = [look_up_table; 0, 0];
elseif strcmp(faces_to_use, 'use_inner_faces')
    mask_1 = all(faces_1 > 0, 2);
    mask_2 = all(faces_2 > 0, 2);
else
    warning('Please specify faces of interest, ''use_complete_faces'' or ''use_all_true_faces'' or ''use_inner_faces'' ');
end

faces_1 = faces_1(mask_1, :);
faces_2 = faces_2(mask_2, :);
faces_1 = sort(faces_1, 2);
faces_2 = sort(faces_2, 2);

% ------------- construct fl_2from1 from fl_1 -------------
map_1to2 = containers.Map(look_up_table(:,1), look_up_table(:,2));
faces_2from1 = - ones(size(faces_1));
for i = 1:size(faces_1, 1)
    for j = 1:2
        if isKey(map_1to2, faces_1(i, j))
            faces_2from1(i, j) = map_1to2(faces_1(i, j));
        end
    end
end

% ------------- sort fl_2from1 and determine a mask by comparing fl_2from1 and fl_2 -------------
faces_2from1 = sort(faces_2from1, 2);
mask = ismember(faces_2from1, faces_2, 'rows');

% ------------- return fl_1(mask, :) and fl_2from1(mask, :) -------------
tracked_uniqueface_1 = faces_1(mask, :);
tracked_uniqueface_2 = faces_2from1(mask, :);

end


