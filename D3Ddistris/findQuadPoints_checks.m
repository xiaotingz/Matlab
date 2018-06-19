% ############################################################################
% This file is several validation checks for findQuadPoints.m, which is
% also related to findTripleLines.m
% ############################################################################

%% ----- 1. #QPs, is it the same as that found by D3D? -----
NType_D3D = double(h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType'))';
NC_D3D = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

NQP_D3D = sum(NType_D3D == 4);
numInnerNode_D3D =  sum(NType_D3D == 2 | NType_D3D == 3 | NType_D3D == 4);
nodeID_D3D = (NC_D3D(:,1)/0.5 + 1) + (NC_D3D(:,2)/0.5)*xDim_N + (NC_D3D(:,3)/0.5)*xDim_N*yDim_N;
QPID_D3D = nodeID_D3D(NType_D3D == 4);
QPID = unique([QPList(:,1);cell2mat(specialQP(:,1))]);

if sum(QPID == QPID_D3D) ~= length(QPID)
    warning('There was some problems of the quad points!');
    QP_both = intersect(QPID, QPID_D3D);
    [QPID_extra, QPID_extra_ind] = setdiff(QPID, QP_both);
    [QPID_D3D_extra, QPID_D3D_extra_ind] = setdiff(QPID_D3D, QP_both);

    % --- check neighboring voxels of the problem node (voxel) from its ID ---
    [I, J, K] = ind2sub([xDim_N, yDim_N, zDim_N], QPID_D3D_extra);
    for idx = 1
        neighbors = [];
        for i = I(idx) - 1 : I(idx)
            for j = J(idx) - 1 : J(idx)
                for k = K(idx) - 1 :K(idx)
                    neighbors = [neighbors, grainIds(i,j,k)];
                end
            end
        end
        unique(neighbors)
    end
end

% ----- 2. How many kinds of specialQPs? What's there physical meaning? How to deal with them? ----- 


% ----- 3. Does the synthetic microstructure satisfy Euler characteristic now? -----






