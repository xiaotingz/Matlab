function [faceCoords, UFlabels] = makeFaceCoords(file)
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

FL = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
    % -- NOTE triNodes are indexes starting from zero 
triNodes = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
NCoords = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

% ##### filter bad data #####
mask = all(FL > 0, 2);
FL = FL(mask,:);
triNodes = triNodes(mask,:);

% ##### sort triangle data (node coords) by facelabel. #####
FL = sort(FL, 2);
data = [FL, zeros(length(FL), 1), triNodes];
data = sortrows(data);

% ##### record facelabels to for tracking of faces #####
UFlabels = [];
% ##### make the {m, n_i, 3} cell #####
idx_Coords = 1;
idx_Info = 1;
faceNodes = [];
faceCoords = {};
for i = 1:length(data)-1
    if data(i, 1) == data(i+1, 1) && data(i, 2) == data(i+1, 2)
        faceNodes = [faceNodes; data(i, 4:6)];
    else
        faceNodes = [faceNodes; data(i, 4:6)];
        faceCoords{idx_Coords} = num2cell(NCoords(faceNodes,:));
        UFlabels = [UFlabels, data(i, 1:2)];
        idx_Coords = idx_Coords + 1;
        faceNodes = [];
    end
end
faceNodes = [faceNodes; data(i, 4:6)];
faceCoords{idx_Coords} = num2cell(NCoords(faceNodes,:));
UFlabels = [UFlabels, data(i, 1:2)];
end










