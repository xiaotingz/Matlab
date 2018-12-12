function [num_nnface_avgcorner] = findGrainFaceNNAvgCorner(file, num_corners)
% ############################################################################
% * Inputs 
%     - num_corners, returned by getFaceCharacter.m
%         Has to be for the full unique faces
% * Notes
%     - Objective: for each grain face, find the average number of corners
%     for each grain face's nearest neighbor faces. Similar to the #neighbor for
%     nearest neighbor of grains.
%     - Dependency: findGrainFaceConnection.m
% ############################################################################
% ----------------------- load debug data -----------------------
% file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
% % '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d'
% faces = faces_an5;
% neigh_list_faceid = neigh_list_faceid_an5;
% num_neigh_face = num_neigh_face_an5;
% num_corners = num_corners_an5;
% ---------------------------------------------------------------

% ##### Find the Connections #####
[faces, num_neigh_face, neigh_list_faceid] = findGrainFaceConnection(file);

% ##### Check Corresp between num_corners and num_neig_face #####
if length(num_corners) ~= length(faces)
    warning('The num_corners input is wrong! Must be for full unique inner grain faces. ')
end

num_nnface_avgcorner = zeros(size(num_corners));
idx = 1;
for i = 1:length(num_corners)
    face_neigh = neigh_list_faceid(idx : idx+num_neigh_face(i)-1);
    nnface_corner = num_corners(face_neigh);
    num_nnface_avgcorner(i) = sum(nnface_corner)/length(nnface_corner);
    idx = idx + num_neigh_face(i);
end

end






