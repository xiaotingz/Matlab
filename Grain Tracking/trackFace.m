function [faces_1_raw, faces_2_raw, corresp] = trackFace(file_1, file_2,look_up_table, use_complete_faces)
% ###########################################################################
% * Input
%     - lookUp = [n,2] = [ID_s1, ID_s2] 
%         is a table which contains only the tracked faces and was sorted in the order of state2.
% * Output
%     - faces = [n,2] = [label_a, label_b] are raw lists including all possible faces. Note just the ones that are tracked.
%     - corresp = [n,2] = [faceID_state1, faceID_state2]
% * Note 
%     - The facesList is duplicated by [A,B] and [B,A]. This code only gives half correspondence. 
%     - This file only asks the grain faces to be complete in their own
%     states but not complete in both states. 
% ###########################################################################
% ------------------ load data for debug --------------------
% file_1 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% file_2 = '/Volumes/XIAOTING/Ni/An5new6_Hsmooth.dream3d';
% use_complete_faces = 0;
% load('look_up_table_an4_an5.mat');
% -----------------------------------------------------------
    num_neigh_1 = double(h5read(file_1,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
    neigh_list_1 = double(h5read(file_1,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
    num_neigh_2 = double(h5read(file_2,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
    neigh_list_2 = double(h5read(file_2,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
    num_neigh_1(1) = [];
    num_neigh_2(1) = [];
    if strcmp(use_complete_faces, 'use_complete_faces')
        surf_grain_1 = h5read(file_1,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
        surf_grain_2 = h5read(file_2,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
        surf_grain_1(1) = [];
        surf_grain_2(1) = [];
    end

    % ##### Find Face Correspondence #####
    %	--- make a look up list that converts ID_An5 to ID_An4. ---
    %	---	if no correspondence, assign NaN ---
    %    *  lookUp_2to1 is made the following way:
    %      -  It was a table of size [N,2], in which N is the number of grains in state_2
    %      -  The first column is the grainID in state_1 and the second column is the grainID in state_2.
    %      -  The rows of the table was sorted so the second column is 1:N. 
    %               Note if a grain had ID_state_2=i but had no correspondence in state_1, the rows is [NaN, i]
    %      -  Now ID_state_2 is implicit. Take the first column to be lookUp_2to1, then it will give the ID_state_1 for grains in state_2
    look_up_table = sortrows(look_up_table, 2);  %  somehow the look up table is not really sorted well 
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

    faces_1_raw = [];
    for i = 1:length(num_neigh_1)
        neighbors = getNeighList(i, num_neigh_1, neigh_list_1);
        faces_1_raw = [faces_1_raw; [i*ones(size(neighbors)), neighbors]];
    end
    faces_1 = [(1:length(faces_1_raw)).', faces_1_raw];
    faces_1 = sortrows(faces_1,[2,3]);
    
    faces_2_raw = [];
    for i = 1:length(num_neigh_2)
        neighbors = getNeighList(i, num_neigh_2, neigh_list_2);
        faces_2_raw = [faces_2_raw; [i*ones(size(neighbors)), neighbors]];
    end
    faces_2_id1 = [(1:length(faces_2_raw))', look_up_table_2to1(faces_2_raw)];
    
    % ##### get rid of the incomplete faces #####
    if strcmp(use_complete_faces, 'use_complete_faces')
        mask_completefaces1 = any(surf_grain_1(faces_1_raw) == 0, 2);
        mask_completefaces2 = any(surf_grain_2(faces_2_raw) == 0, 2);
        faces_1 = faces_1(mask_completefaces1, :);
        faces_2_id1 = faces_2_id1(mask_completefaces2, :);
    end
    
    % ##### sort faces_state2 into ID of state1 #####
    mask = ~(any(isnan(faces_2_id1),2));
    faces_2_id1 = faces_2_id1(mask,:);
    faces_2_id1 = sortrows(faces_2_id1,[2,3]);

    % ##### compare the common element thus track the faces, when FULL list #####
    length_1 = length(faces_1);
    length_2 = length(faces_2_id1);
    i = 1;
    j = 1;
    corresp = [];
    while j <= length_1 && i <= length_2
        if faces_2_id1(i,2) == faces_1(j,2)
            if faces_2_id1(i,3) == faces_1(j,3)
                corresp = vertcat(corresp, [faces_1(j,1), faces_2_id1(i,1)]);
                i = i + 1;
                j = j + 1;
            elseif faces_2_id1(i,3) > faces_1(j,3)
                j = j + 1;
            else
                i = i + 1;
            end
        elseif faces_2_id1(i,2) > faces_1(j,2)
            j = j + 1;
        else
            i = i + 1;
        end
    end
    corresp = sortrows(corresp, 2);
end
