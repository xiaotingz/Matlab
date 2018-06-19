function result = findQuadPoints(file)
% ############################################################################
% - QPs = [n, 4] = [G1, G2, G3, G4]
%     G1-G4 are the four grains encapsulating the point
% - QPpos = [n, 3] = [x, y, z]
%     x, y, z are the QuadPoint positions, max(x) = numVoxel_x + 1
% - The algorithm is inspired from D3D's quick mesh.
%     - I care only the inner volume QuadPoints. 
%     - Every voxel has 6 neighbors, or 6 faces, each face have 4 nodes.
%     To avoid double counting, only need to check 3 neighbors of each
%     voxel, except for the initial ones.
%     - Also, see note Topology written by GoodNotes
% ############################################################################
% ----------------------- load debug data -----------------------
% clear
% file = '/Users/xiaotingzhong/Desktop/Datas/synthetic/180502_CubicSingleEquiaxedOut.dream3d';
% ---------------------------------------------------------------
dims = double(h5read(file, '/DataContainers/SyntheticVolumeDataContainer/_SIMPL_GEOMETRY/DIMENSIONS'));
grainIds = squeeze(h5read(file, '/DataContainers/SyntheticVolumeDataContainer/CellData/FeatureIds'));

% --- dims are the dimension of voxels, dimension of nodes is (dims + 1) ---
xDim_N = dims(1) + 1;
yDim_N = dims(2) + 1;
zDim_N = dims(3) + 1;
layerSize_N = xDim_N * yDim_N;
ownerList = cell(xDim_N*yDim_N*zDim_N, 1);
for i = 2 : (dims(1) - 1)
    for j = 2 : (dims(2) - 1)
        for k = 2 : (dims(3)- 1)

            voxel = grainIds(i, j, k);
            neigh1 = grainIds(i+1, j, k);
            neigh2 = grainIds(i, j+1, k);
            neigh3 = grainIds(i, j, k+1);
            
            if i == 2
                neigh4 = grainIds(i-1, j, k);
                if voxel ~= neigh4
                    node1 =  calcNodeID(i, j, k, xDim_N, yDim_N);
                    node2 =  calcNodeID(i, j, k+1, xDim_N, yDim_N);
                    node3 =  calcNodeID(i, j+1, k, xDim_N, yDim_N);
                    node4 =  calcNodeID(i, j+1, k+1, xDim_N, yDim_N);
                    ownerList{node1} = insert(ownerList{node1}, [voxel, neigh4]);    
                    ownerList{node2} = insert(ownerList{node2}, [voxel, neigh4]);
                    ownerList{node3} = insert(ownerList{node3}, [voxel, neigh4]);            
                    ownerList{node4} = insert(ownerList{node4}, [voxel, neigh4]);          
                end
            end
            if j == 2
                neigh5 = grainIds(i, j-1, k);
                if voxel ~= neigh5
                    node1 =  calcNodeID(i+1, j, k, xDim_N, yDim_N);
                    node2 =  calcNodeID(i+1, j, k+1, xDim_N, yDim_N);
                    node3 =  calcNodeID(i, j, k, xDim_N, yDim_N);
                    node4 =  calcNodeID(i, j, k+1, xDim_N, yDim_N);
                    ownerList{node1} = insert(ownerList{node1}, [voxel, neigh5]);    
                    ownerList{node2} = insert(ownerList{node2}, [voxel, neigh5]);
                    ownerList{node3} = insert(ownerList{node3}, [voxel, neigh5]);            
                    ownerList{node4} = insert(ownerList{node4}, [voxel, neigh5]);          
                end
            end
            if k == 2
                neigh6 = grainIds(i, j, k-1);
                if voxel ~= neigh6
                    node1 =  calcNodeID(i+1, j, k, xDim_N, yDim_N);
                    node2 =  calcNodeID(i+1, j+1, k, xDim_N, yDim_N);
                    node3 =  calcNodeID(i, j, k, xDim_N, yDim_N);
                    node4 =  calcNodeID(i, j+1, k, xDim_N, yDim_N);
                    ownerList{node1} = insert(ownerList{node1}, [voxel, neigh6]);    
                    ownerList{node2} = insert(ownerList{node2}, [voxel, neigh6]);
                    ownerList{node3} = insert(ownerList{node3}, [voxel, neigh6]);            
                    ownerList{node4} = insert(ownerList{node4}, [voxel, neigh6]);          
                end
            end
            
            if voxel ~= neigh1
                node1 =  calcNodeID(i+1, j, k, xDim_N, yDim_N);
                node2 =  calcNodeID(i+1, j, k+1, xDim_N, yDim_N);
                node3 =  calcNodeID(i+1, j+1, k, xDim_N, yDim_N);
                node4 =  calcNodeID(i+1, j+1, k+1, xDim_N, yDim_N);
                ownerList{node1} = insert(ownerList{node1}, [voxel, neigh1]);    
                ownerList{node2} = insert(ownerList{node2}, [voxel, neigh1]);
                ownerList{node3} = insert(ownerList{node3}, [voxel, neigh1]);            
                ownerList{node4} = insert(ownerList{node4}, [voxel, neigh1]);          
            end
            if voxel ~= neigh2
                node1 =  calcNodeID(i, j+1, k, xDim_N, yDim_N);
                node2 =  calcNodeID(i, j+1, k+1, xDim_N, yDim_N);
                node3 =  calcNodeID(i+1, j+1, k, xDim_N, yDim_N);
                node4 =  calcNodeID(i+1, j+1, k+1, xDim_N, yDim_N);
                ownerList{node1} = insert(ownerList{node1}, [voxel, neigh2]);    
                ownerList{node2} = insert(ownerList{node2}, [voxel, neigh2]);
                ownerList{node3} = insert(ownerList{node3}, [voxel, neigh2]);            
                ownerList{node4} = insert(ownerList{node4}, [voxel, neigh2]);          
            end
            if voxel ~= neigh3
                node1 =  calcNodeID(i, j, k+1, xDim_N, yDim_N);
                node2 =  calcNodeID(i, j+1, k+1, xDim_N, yDim_N);
                node3 =  calcNodeID(i+1, j, k+1, xDim_N, yDim_N);
                node4 =  calcNodeID(i+1, j+1, k+1, xDim_N, yDim_N);
                ownerList{node1} = insert(ownerList{node1}, [voxel, neigh3]);    
                ownerList{node2} = insert(ownerList{node2}, [voxel, neigh3]);
                ownerList{node3} = insert(ownerList{node3}, [voxel, neigh3]);            
                ownerList{node4} = insert(ownerList{node4}, [voxel, neigh3]);          
            end
        end
    end
end

% ##### Make lists for the Quads & superQuads #####
% ----- Quads -----
QPList = [];
superQPs = {};
idx_QP = 1;
idx_SQP = 1;
for i = 1:length(ownerList)
    if ~isempty(ownerList{i})
        if size(ownerList{i},2) == 4
            QPList(idx_QP, 1) = i;
            QPList(idx_QP, 2:5) = cell2mat(ownerList{i});
            idx_QP = idx_QP + 1;
        elseif size(ownerList{i},2) > 4
            superQPs{idx_SQP, 1} = i;
            superQPs{idx_SQP, 2} = ownerList{i};
            idx_SQP = idx_SQP + 1;
        end
    end
end

% ----- super Quads -----
FivePList = [];
SixPList = [];
SevenPList = [];
EightPList = [];
idx_FC = 1;
idx_SiC = 1;
idx_SeC = 1;
idx_EC = 1;

for i = 1:length(superQPs)
    if length(superQPs{i,2}) == 5
        FivePList(idx_FC, 1) = superQPs{i,1};
        FivePList(idx_FC, 2:6) = cell2mat(superQPs{i,2});
        idx_FC = idx_FC + 1;
    else
        if length(superQPs{i,2}) == 6
            SixPList(idx_SiC, 1) = superQPs{i,1};
            SixPList(idx_SiC, 2:7) = cell2mat(superQPs{i,2});
            idx_SiC = idx_SiC + 1;
        elseif length(superQPs{i,2}) == 7
            SevenPList(idx_SeC, 1) = superQPs{i,1};
            SevenPList(idx_SeC, 2:8) = cell2mat(superQPs{i,2});
            idx_SeC = idx_SeC + 1;
        elseif length(superQPs{i,2}) == 8
            EightPList(idx_EC, 1) = superQPs{i,1};
            EightPList(idx_EC, 2:9) = cell2mat(superQPs{i,2});
            idx_EC = idx_EC + 1;
        end
    end
end
result = {QPList, FivePList, SixPList, SevenPList, EightPList};

% ##### Get coordinates from nodeID/voxelID #####
% [I, J, K] = ind2sub([xDim_N, yDim_N, zDim_N], nodeID);

end






