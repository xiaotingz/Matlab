function [num_nnface_avgcorner] = findFaceNNAvgCorner(file, faces, num_corners, tl)
% ############################################################################
% * Inputs 
%     - num_corners [n, 1] & faces [n, 2], returned by getFaceCharacter.m
%         MUST BE for the FULL unique faces
%     - tl = [m, 3]
%         Returned by de findTripleLines.m. Also, don't threshold triangle
%         quality for connection purpose
% * Notes
%     - Objective: for each grain face, find the average number of corners
%     for each grain face's nearest neighbor faces. Similar to the #neighbor for
%     nearest neighbor of grains.
%     - Dependency: findFaceConnection.m
% ############################################################################
% ----------------------- load debug data -----------------------
% file = file_an4;
% tl = triple_line_full_an4;
% faces = faces_an4;
% num_corners = num_corners_an4;
% ---------------------------------------------------------------

% ##### Find the Connections #####
% """
% Note this is for unique faces, not the full faces.
% """
[faces_2, num_neigh_face, neigh_list_faceid] = findFaceConnection(file, tl);

% -----------------------------------------------------------------------------
% % % % ##### OLD VERSION, Check Corresp between num_corners and num_neig_face #####
% % % if any(sum(faces_2 ~= faces) > 0, 2)
% % %     warning('The num_corners input is wrong! Must be for full unique inner grain faces. ')
% % % end
% % 
% % % num_nnface_avgcorner = zeros(size(num_corners));
% % % idx = 1;
% % % for i = 1:length(num_corners)
% % %     face_neigh = neigh_list_faceid(idx : idx+num_neigh_face(i)-1);
% % %     nnface_corner = num_corners(face_neigh);
% % %     num_nnface_avgcorner(i) = sum(nnface_corner)/length(nnface_corner);
% % %     idx = idx + num_neigh_face(i);
% % % end
% -----------------------------------------------------------------------------

% ##### Find Connections For The Full Faces #####
% ----- build a dictionary from the unique faces -----
keys_faces = num2str(faces_2);
keys_faces = mat2cell(keys_faces, ones(1, size(keys_faces, 1)), size(keys_faces, 2));
faceid_dict = containers.Map(keys_faces, (1:length(faces_2))');

% ----- convert num_corners to a unique list -----
% """
% This step is redundant if num_corners was found for the unique faces,
% necessary if full faces are passed in.
% """
faces = sort(faces, 2);
[~, idx_unique, ~] = unique(faces, 'rows');
tmp = faces(idx_unique, :);
if any(sum(faces_2 ~= tmp) > 0, 2)
    warning('The corresp between num_corners & face connections are WRONG! ')
end
num_corners_uniqueface = num_corners(idx_unique);

% ----- build a set of keys from full faces, which is consistent with num_coners -----
num_nnface_avgcorner = zeros(size(faces, 1), 1);
faces_keys = num2str(faces);
for i = 1:length(faces_keys)
    face_id = faceid_dict(faces_keys(i, :));
    face_neigh = getNeighList(face_id, num_neigh_face, neigh_list_faceid);
    nnface_corner = num_corners_uniqueface(face_neigh);
    num_nnface_avgcorner(i) = sum(nnface_corner)/length(nnface_corner);
end

end






