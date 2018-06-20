function result = findQuadPoints(file)
% ############################################################################
% - QNs = [n, 4] = [G1, G2, G3, G4]
%     G1-G4 are the four grains encapsulating the point
% - QNpos = [n, 3] = [x, y, z]
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
% file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_mesh.dream3d');
% ---------------------------------------------------------------
% dims = double(h5read(file, '/DataContainers/SyntheticVolumeDataContainer/_SIMPL_GEOMETRY/DIMENSIONS'));
% grainIds = squeeze(h5read(file, '/DataContainers/SyntheticVolumeDataContainer/CellData/FeatureIds'));
dims = double(h5read(file, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS'));
grainIds = squeeze(h5read(file, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));
surfG = boolean(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures'))';
surfG(1) = [];

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
QNList = [];
superQNs = {};
idx_QN = 1;
idx_SQN = 1;
for i = 1:length(ownerList)
    if ~isempty(ownerList{i})
        if size(ownerList{i},2) == 4
            QNList(idx_QN, 1) = i;
            QNList(idx_QN, 2:5) = cell2mat(ownerList{i});
            idx_QN = idx_QN + 1;
        elseif size(ownerList{i},2) > 4
            superQNs{idx_SQN, 1} = i;
            superQNs{idx_SQN, 2} = ownerList{i};
            idx_SQN = idx_SQN + 1;
        end
    end
end

% ----- super Quads -----
fiveCoordNList = [];
sixCoordNList = [];
sevenCoordNList = [];
eightCoordNList = [];
idx_FC = 1;
idx_SiC = 1;
idx_SeC = 1;
idx_EC = 1;

for i = 1:length(superQNs)
    if length(superQNs{i,2}) == 5
        fiveCoordNList(idx_FC, 1) = superQNs{i,1};
        fiveCoordNList(idx_FC, 2:6) = cell2mat(superQNs{i,2});
        idx_FC = idx_FC + 1;
    else
        if length(superQNs{i,2}) == 6
            sixCoordNList(idx_SiC, 1) = superQNs{i,1};
            sixCoordNList(idx_SiC, 2:7) = cell2mat(superQNs{i,2});
            idx_SiC = idx_SiC + 1;
        elseif length(superQNs{i,2}) == 7
            sevenCoordNList(idx_SeC, 1) = superQNs{i,1};
            sevenCoordNList(idx_SeC, 2:8) = cell2mat(superQNs{i,2});
            idx_SeC = idx_SeC + 1;
        elseif length(superQNs{i,2}) == 8
            eightCoordNList(idx_EC, 1) = superQNs{i,1};
            eightCoordNList(idx_EC, 2:9) = cell2mat(superQNs{i,2});
            idx_EC = idx_EC + 1;
        end
    end
end

% ##### for pillar shape sample, starting with i,j,k=2 won't exclude the surface voxels #####
QNs = unique(QNList(:,2:5), 'rows');
% --- if one of the ID=0, the voxel is on freeSurface, thus doesn't give a quad --- 
QNs = QNs(~any(QNs == 0, 2), :);
% --- if all the 4 grains are surface grains, the voxel doesn't give a quad --- 
QNs = QNs(sum(surfG(uint8(QNs)),2) ~= 4, :);

% --- same process applies to the superQuads, note 0 can only show once in the first digit --- 
fiveCoordNs = []; sixCoordNs = []; sevenCoordNs = []; eightCoordNs = [];
if ~isempty(fiveCoordNList)
    fiveCoordNs = unique(fiveCoordNList(:,2:6), 'rows');
    mask_QNinFN = any(fiveCoordNs == 0, 2);
    QNs = [QNs; fiveCoordNs(mask_QNinFN, 2:5)];
    fiveCoordNs = fiveCoordNs(~mask_QNinFN, :);
    fiveCoordNs = fiveCoordNs(sum(surfG(uint8(fiveCoordNs)),2) ~= 5, :);
end
if ~isempty(sixCoordNList)
    sixCoordNs = unique(sixCoordNList(:,2:7), 'rows');
    mask_FNinSixN = any(sixCoordNs == 0, 2);
    fiveCoordNs = [fiveCoordNs; sixCoordNs(mask_FNinSixN, 2:6)];
    sixCoordNs = sixCoordNs(~mask_FNinSixN, :);
    sixCoordNs = sixCoordNs(sum(surfG(uint8(sixCoordNs)),2) ~= 6, :);
end
if ~isempty(sevenCoordNList)
    sevenCoordNs = unique(sevenCoordNList(:,2:8), 'rows');
    mask_sixNinSevenN = any(sevenCoordNs == 0, 2);
    sixCoordNs = [sixCoordNs; sevenCoordNs(mask_sixNinSevenN, 2:7)];
    sevenCoordNs = sevenCoordNs(~mask_sixNinSevenN, :);
    sevenCoordNs = sevenCoordNs(sum(surfG(uint8(sevenCoordNs)),2) ~= 7, :);
end
if ~isempty(eightCoordNList)
    eightCoordNs = unique(eightCoordNList(:,2:9), 'rows');
    mask_sevenNinEN = any(eightCoordNs == 0, 2);
    sevenCoordNs = [sevenCoordNs; eightCoordNs(mask_sevenNinEN, 2:8)];
    eightCoordNs = eightCoordNs(~mask_sevenNinEN, :);
    eightCoordNs = eightCoordNs(sum(surfG(uint8(eightCoordNs)),2) ~= 8, :);
end
result = {QNs, fiveCoordNs, sixCoordNs, sevenCoordNs, eightCoordNs};

% QNs = unique(QNList(:,2:5), 'rows');
% 
% fiveCoordNs = []; sixCoordNs = []; sevenCoordNs = []; eightCoordNs = [];
% if ~isempty(FiveNList)
%     fiveCoordNs = unique(FiveNList(:,2:6), 'rows');
% end
% if ~isempty(SixNList)
% sixCoordNs = unique(SixNList(:,2:7), 'rows');
% end
% if ~isempty(SevenNList)
% sevenCoordNs = unique(SevenNList(:,2:8), 'rows');
% end
% if ~isempty(EightNList)
% eightCoordNs = unique(EightNList(:,2:9), 'rows');
% end
% result = {QNs, fiveCoordNs, sixCoordNs, sevenCoordNs, eightCoordNs};

% ##### Get coordinates from nodeID/voxelID #####
% [I, J, K] = ind2sub([xDim_N, yDim_N, zDim_N], nodeID);

end





