function [face_coords, unique_facelabel] = makeFaceCoords(file)
% ##########################################################################
% * Input
%     - faces([N,2]) and faceCorresp([N,1]) is returned by the function 'TrackFace'
% * Output
%     - FaceCoords = {m, n_i, 3} 
%         m: #uniqueFaces. 
%               Note this is different from all other function in which
%               [A, B] and [B, A] are copied to be the same because cell is slow. 
%         n_i:  #nodes on the ith face
%         3: x, y, z of each node 
%     - FaceInfo = [N, 3] = [idx, FL_A, FL_B]
%         idx, is the index that converst facesCoords into the order of faces. 
%         The procedure is get faceM(matrix of [N,3]) from faceCoords, then 
%             faceM = [faceM, faceM]';
%             reshape(faceM, [], 1);
%             faceM = faceM(faceInfo(:,1));
% * NOTE
%     - make sure samples are aligned by changing ORIGIN.
% ##########################################################################
% ----------------------- load debug data -----------------------
% file = file_An5;
% faces = faces_An5;
% ---------------------------------------------------------------

facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
    % -- NOTE triNodes are indexes starting from zero 
tri_node = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

% ##### filter bad data #####
mask = all(facelabel > 0, 2);
facelabel = facelabel(mask,:);
tri_node = tri_node(mask,:);

% ##### sort triangle data (node coords) by facelabel. #####
facelabel = sort(facelabel, 2);
data = [facelabel, zeros(length(facelabel), 1), tri_node];
data = sortrows(data);

% ##### record facelabels to for tracking of faces #####
unique_facelabel = [];
% ##### make the {m, n_i, 3} cell #####
idx_coords = 1;
face_nodes = [];
face_coords = {};
for i = 1:length(data)-1
    if data(i, 1) == data(i+1, 1) && data(i, 2) == data(i+1, 2)
        face_nodes = [face_nodes; data(i, 4:6)];
    else
        face_nodes = [face_nodes; data(i, 4:6)];
        face_coords{idx_coords} = num2cell(node_coord(face_nodes,:));
        unique_facelabel = [unique_facelabel, data(i, 1:2)];
        idx_coords = idx_coords + 1;
        face_nodes = [];
    end
end
face_nodes = [face_nodes; data(i, 4:6)];
face_coords{idx_coords} = num2cell(node_coord(face_nodes,:));
unique_facelabel = [unique_facelabel, data(i, 1:2)];
end










