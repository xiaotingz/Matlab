function nodeID = calcNodeID(i,j,k, xDim_N, yDim_N)
% ############################################################################
% This is an auxiliary function for findQuadPints
% - Node xDim_N is the dimension of the nodes, which equals xDim_sample + 1
% ############################################################################
    nodeID = i + (j-1)*xDim_N + (k-1)*xDim_N*yDim_N;
end