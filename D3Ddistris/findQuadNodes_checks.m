% ############################################################################
% This file is several validation checks for findQuadPoints.m, which is
% also related to findTripleLines.m
%  1. #QNs, is it the same as that found by D3D?
%  2. Does the synthetic microstructure satisfy Euler characteristic now?
% ############################################################################

% %% ###### 1. #QNs, is it the same as that found by D3D? ######
file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_mesh.dream3d');
dims = double(h5read(file, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS'));
xDim_N = dims(1) + 1;
yDim_N = dims(2) + 1;
zDim_N = dims(3) + 1;

NType_D3D = double(h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType'))';
NC_D3D = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

NQN_D3D = sum(NType_D3D == 4);
numInnerNode_D3D =  sum(NType_D3D == 2 | NType_D3D == 3 | NType_D3D == 4);
nodeID_D3D = (NC_D3D(:,1)/0.5 + 1) + (NC_D3D(:,2)/0.5)*xDim_N + (NC_D3D(:,3)/0.5)*xDim_N*yDim_N;
QNID_D3D = nodeID_D3D(NType_D3D == 4);
QNID = unique([QNList(:,1); fiveCoordNList(:,1); sixCoordNList(:,1)]);

if sum(QNID == QNID_D3D) ~= length(QNID)
    warning('There was some problems of the quad points!');
    QN_both = intersect(QNID, QNID_D3D);
    [QNID_extra, QNID_extra_ind] = setdiff(QNID, QN_both);
    [QNID_D3D_extra, QNID_D3D_extra_ind] = setdiff(QNID_D3D, QN_both);

    % --- check neighboring voxels of the problem node (voxel) from its ID ---
    [I, J, K] = ind2sub([xDim_N, yDim_N, zDim_N], QNID_D3D_extra);
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

% ###### 2. How many kinds of specialQNs? What's there physical meaning? How to deal with them? ######


%% ###### 3. Does the synthetic microstructure satisfy Euler characteristic now? ######
file = '/Users/xiaotingzhong/Desktop/Datas/synthetic/180502_CubicSingleEquiaxedOut.dream3d';

numNeigh = double(h5read(file,'/DataContainers/SyntheticVolumeDataContainer/CellFeatureData/NumNeighbors'))';
numNeigh(1) = [];
NeighborList = double(h5read(file,'/DataContainers/SyntheticVolumeDataContainer/CellFeatureData/NeighborList'));
surfG = boolean(h5read(file,'/DataContainers/SyntheticVolumeDataContainer/CellFeatureData/SurfaceFeatures'))';
surfG(1) = [];


F = numNeigh;

E = zeros(size(F));
V_SQasNormal = zeros(size(F));
V_SQNlusOne = zeros(size(F));

% QNs = result{1,1};
% FCNs = result{1,2};

for i = 1:length(TLs)
    for j = 1:3
        E(TLs(i,j)) = E(TLs(i,j)) + 1;
    end
end

for i = 1:length(QNs)
    for j = 1:4
        V_SQasNormal(QNs(i,j)) = V_SQasNormal(QNs(i,j)) + 1;
        V_SQNlusOne(QNs(i,j)) = V_SQNlusOne(QNs(i,j)) + 1;
    end
end
for i = 1:length(FCNs)
    for j = 1:5
        V_SQasNormal(FCNs(i,j)) = V_SQasNormal(FCNs(i,j)) + 1;
        V_SQNlusOne(FCNs(i,j)) = V_SQNlusOne(FCNs(i,j)) + 2;
    end
end

eulerChar = [(1:length(F))', F, E, V_SQasNormal, V_SQasNormal-E+F];
eulerChar = eulerChar(~surfG, :);
generalQNs = [[QNs, zeros(length(QNs),1)]; FCNs];

%%
objGrain = 16;
neighbors = getNeighList(objGrain, numNeigh, NeighborList)
% --- generalQNs = [n, c], see function getFaceCharacter for details --- 

for i = 1:length(neighbors)
    [numCorners, numEdges] = getFaceCharacter([objGrain, neighbors(i)], TLs, generalQNs);
    disp(['for [', num2str(objGrain), ', ', num2str(neighbors(i)), ']: numCorners = ', num2str(numCorners), ', numEdges = ', num2str(numEdges)]);
end

%% --- check the TLs and Quads of a specific grainFace ---
idx = 10;
specialFace(idx, :)
A = specialFace(idx, 1);
B = specialFace(idx, 2);
mask_TL_A = (TLs == A);
mask_TL_B = (TLs == B);
mask_TL_AB = ((sum(mask_TL_A, 2) + sum(mask_TL_B, 2)) == 2);
TLs(mask_TL_AB, :);
mask_QN_A = (generalQNs == A);
mask_QN_B = (generalQNs == B);
mask_QN_AB = ((sum(mask_QN_A, 2) + sum(mask_QN_B, 2)) == 2);
generalQNs(mask_QN_AB, :)


%% ##### superQuad multiplicity check #####
% --- If a FiveCoordNode always corresponds to two quads, 
% --- then to get the correct #Quads, one can simply ignore the FiveCoordNodes
% labels = [zeros(length(NeighborList),1), NeighborList];
% cnt = 0;
% for i = 1:length(numNeigh)
%     labels(cnt+1 : cnt+numNeigh(i), 1) = ones(numNeigh(i),1)*i;
%     cnt = cnt + numNeigh(i);
% end

innerGs = [(1:length(numNeigh))', numNeigh];
innerGs = innerGs(~surfG,:);
labels = zeros(sum(innerGs(:,2)),2);
cnt = 0;
for i = 1:length(innerGs)
    labels(cnt+1 : cnt+innerGs(i,2), 1) = ones(innerGs(i,2),1) * innerGs(i,1);
    labels(cnt+1 : cnt+innerGs(i,2), 2) = getNeighList(innerGs(i,1), numNeigh, NeighborList);
    cnt = cnt + innerGs(i,2);
end

[numCorners, numEdges] = getFaceCharacter(labels, TLs, QNs, FCNs);

mask_CeqE = (numCorners ~= numEdges);
unbalanceFace = [labels, numCorners, numEdges];
unbalanceFace = unbalanceFace(mask_CeqE, :);

mask_specialFace = (unbalanceFace(:,3) - unbalanceFace(:,4) ~= 1);
specialFace = unbalanceFace(mask_specialFace, :);
%%
% ##### clear multiplicityQNs #####
for i = 1:length(FNmulti)
    mask_multi = (FCNs(i,1)==QNs |FCNs(i,2)==QNs |FCNs(i,3)==QNs |FCNs(i,4)==QNs |FCNs(i,5)==QNs);
    mask_multi = (sum(mask_multi,2) == 4);
    QNs(mask_multi,:) = [];
end



% ###### 3. How many kinds of specialQNs? What's there physical meaning? How to deal with them? ######

