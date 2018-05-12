% function numEdges = calcFaceNumEdges(file, faces)
function numEdges = calcFaceNumEdges(FL, triNodes, Ntypes, faces)
% ############################################################################
% NOTES
% - assuming the non-function 0 at the first element of numNeighbors has been deleted.
% ############################################################################
% numNeigh(1) = [];
% ----------------------- load debug data -----------------------
% faces = faces_An4;
% file = file_An4;
% ---------------------------------------------------------------

% FL = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
%     % -- NOTE triNodes are indexes starting from zero 
% triNodes = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
% Ntypes = double(h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType'))';
mask_surf = any(FL <=0, 2);
FL(mask_surf, :) = [];
triNodes(mask_surf, :) = [];

numEdges = zeros(length(faces), 1);
for i = 1:length(faces)
    GA = faces(i,1);
    GB = faces(i,2);
    mask_face = ((FL(:,1) == GA & FL(:,2) ==GB ) | (FL(:,1) == GB & FL(:,2) ==GA ));
    faceNodes = triNodes(mask_face, :);
    mask_quads = (Ntypes(faceNodes) == 4);
    faceQuads = unique(faceNodes(mask_quads));
    % ----- always use GA as the grain of interest. !!! MAY BE WRONG !!! -----
    j = 1;
    while j <= length(faceQuads)
        mask_thisQuad = any(triNodes == faceQuads(j), 2);
        labels = FL(mask_thisQuad, :);
        mask_centerGF = any(labels == GA, 2);
        tmp = labels(mask_centerGF,:);
        if length(unique(tmp)) == 4
            j = j +1;
        else
            faceQuads(j) = [];
        end
    end
    numEdges(i) = length(faceQuads); 
end

end




