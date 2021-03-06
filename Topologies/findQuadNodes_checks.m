% ############################################################################
% This file is several validation checks for findQuadPoints.m, which is
% also related to findTripleLines.m
%  1. #QNs, is it the same as that found by D3D?
%  2. Does the synthetic microstructure satisfy Euler characteristic now?
%     - Grain balance check
%     - Grain face balance check: faces of one grain & general faces
%  3. What info should be used to define a good QNs?
%     - #Grains at the QN == 4?
%     - #TLs at the QN == 4?
%     - #GBs at the QN == 6?
%     - How many kinds of specialQNs? What's there physical meaning? How to deal with them? 
% ############################################################################
%% --------------------------- Data preparation ----------------------------
% load('180625_Ni_TopologyResult.mat');
% file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_mesh_mergeTwin.dream3d';
% fiveCoordNList = fiveCoordNList_An4;
% QNList = QNList_An4;
% result = result_An4;
% sixCoordNList = sixCoordNList_An4;
% TLs = TLs_An4;
% % 
load('180625_Ni_TopologyResult.mat');
file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh_mergeTwin.dream3d';
fiveCoordNList = fiveCoordNList_An5;
QNList = QNList_An5;
result = result_An5;
sixCoordNList = sixCoordNList_An5;
TLs = TLs_An5;

clearvars -except file fiveCoordNList QNList result sixCoordNList TLs
%  --------------------------------------------------------------------------


%% ###### 1. #QNs, is it the same as that found by D3D? ######
% file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_mesh.dream3d');
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



%% ###### 2.1. Grain Balance Check: Euler Characteristic ######
% numNeigh = double(h5read(file,'/DataContainers/SyntheticVolumeDataContainer/CellFeatureData/NumNeighbors'))';
% numNeigh(1) = [];
% NeighborList = double(h5read(file,'/DataContainers/SyntheticVolumeDataContainer/CellFeatureData/NeighborList'));
% surfG = boolean(h5read(file,'/DataContainers/SyntheticVolumeDataContainer/CellFeatureData/SurfaceFeatures'))';
% surfG(1) = [];
numNeigh = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors'))';
numNeigh(1) = [];
NeighborList = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
surfG = boolean(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures'))';
surfG(1) = [];


F = numNeigh;
E = zeros(size(F));
V_SQasNormal = zeros(size(F));
V_SQNplusOne = zeros(size(F));
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
        V_SQNplusOne(QNs(i,j)) = V_SQNplusOne(QNs(i,j)) + 1;
    end
end
for i = 1:length(FCNs)
    for j = 1:5
        V_SQasNormal(FCNs(i,j)) = V_SQasNormal(FCNs(i,j)) + 1;
        V_SQNplusOne(FCNs(i,j)) = V_SQNplusOne(FCNs(i,j)) + 2;
    end
end

eulerChar_FCNplus1 = [(1:length(F))', F, E, V_SQNplusOne, V_SQNplusOne-E+F];
eulerChar_FCNplus1 = eulerChar_FCNplus1(~surfG, :);


% % ----- For one grain: balance of its faces -----
% objGrain = 16;
% neighbors = getNeighList(objGrain, numNeigh, NeighborList);
% 
% % """
% % getFaceCharacter(labels, TLs, QNs, FCNs). The FCNs input can be empty.
% % """
% for i = 1:length(neighbors)
%     [numCorners, numEdges] = getFaceCharacter([objGrain, neighbors(i)], TLs, QNs, FCNs);
%     disp(['for [', num2str(objGrain), ', ', num2str(neighbors(i)), ']: numCorners = ', num2str(numCorners), ', numEdges = ', num2str(numEdges)]);
% end

%% ##### 2.1. Plot ##### 
% ----- EC = [(1:length(F))', F, E, V, V-E+F] ----- 
histogram(EC_An4(:,5),'Normalization','probability')
hold on
histogram(EC_An5(:,5),'Normalization','probability')
xlim([-5, 15])
legend('An4','An5','Location','NW');

xlabel('V - E + F');
ylabel('Count');
title('using cleanedGoodQuads only')

print('EulerCharacteristics', '-dpng','-r300')


%% ##### 2.2. Face Balance Check: Corner & Edge  #####
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

%  ----- Fone one grain face, get the TLs and Quads -----
idx = 10;
specialFace(idx, :)
A = specialFace(idx, 1);
B = specialFace(idx, 2);
mask_TL_A = (TLs == A);
mask_TL_B = (TLs == B);
mask_TL_AB = ((sum(mask_TL_A, 2) + sum(mask_TL_B, 2)) == 2);
TLs(mask_TL_AB, :);
mask_QN_A = (QNs == A);
mask_QN_B = (QNs == B);
mask_QN_AB = ((sum(mask_QN_A, 2) + sum(mask_QN_B, 2)) == 2);
QNs(mask_QN_AB, :)

%% ##### 3.1. Remove the multiplicity QNs from FCNs, then find #TLs and #GBs at the nodes #####
% """
% - Each FCN defines a five-grain set, which probably stemed from insufficient experimental resolution.
% - Then within each five-grain set, there should be 2 and only 2 QNs. But that's not the case, 
%       probably because of artifact from reconstruction & voxelizaiton
% - Decided to first remove all QNs stemmed from FCNs (muliplicities). 
% - Then do calculations from the non-multiple-QNs & FCNs. 
% """

% FCNs = fiveCoordNList;

QNs = result{1,1};
FCNs = result{1,2};
FCN_multi = zeros(length(FCNs),1);
for i = 1:length(FCNs)
    mask_multi = (FCNs(i,1)==QNs |FCNs(i,2)==QNs |FCNs(i,3)==QNs |FCNs(i,4)==QNs |FCNs(i,5)==QNs);
    mask_multi = (sum(mask_multi,2) == 4);
    FCN_multi(i) = sum(mask_multi);
    QNs(mask_multi,:) = [];
end

SixCNs = result{1,3};
SixCN_multiFour = zeros(length(SixCNs),1);
for i = 1:length(SixCNs)
    mask_multiFour = (SixCNs(i,1)==QNs |SixCNs(i,2)==QNs |SixCNs(i,3)==QNs |SixCNs(i,4)==QNs |SixCNs(i,5)==QNs |SixCNs(i,6)==QNs);
    mask_multiFour = (sum(mask_multiFour,2) == 4);
    SixCN_multiFour(i) = sum(mask_multiFour);
    QNs(mask_multiFour,:) = [];
end

%% ------ #TLs and #GBs -----
% """
% objType = 'QNs' | 'FCNs' | 'SixCNs'
% """
SixCNs = result{1,3};

[numGBs_QN, numTLs_QN] = nodeInfo('QNs', file, TLs, QNs, FCNs, SixCNs);

fewGB_QN = [(1:length(numGBs_QN))', numGBs_QN];
fewGB_QN = fewGB_QN(numGBs_QN~=6, :);
fewTL_QN = [(1:length(numTLs_QN))', numTLs_QN];
fewTL_QN = fewTL_QN(numTLs_QN~=4, :);

QNs_An5= QNs;
fewGB_QN_An5 = fewGB_QN;
fewTL_QN_An5 = fewTL_QN;
save('An5_QN.mat', 'QNs_An5','fewGB_QN_An5','fewTL_QN_An5');

% [numGBs_FCN, numTLs_FCN] = nodeInfo('FCNs', file, TLs, QNs, FCNs, SixCNs);
% 
% [numGBs_SixCN, numTLs_SixCN] = nodeInfo('SixCNs', file, TLs, QNs, FCNs, SixCNs);


%% ##### 3.2.  SixCN multiplicity from non-cleared QNs #####
QNs = result{1,1};
FCNs = result{1,2};
SixCNs = result{1,3};
SixCN_multiFour = zeros(length(SixCNs),1);
SixCN_multiFive = zeros(length(SixCNs),1);
for i = 1:length(SixCNs)
    mask_multiFour = (SixCNs(i,1)==QNs |SixCNs(i,2)==QNs |SixCNs(i,3)==QNs |SixCNs(i,4)==QNs |SixCNs(i,5)==QNs |SixCNs(i,6)==QNs);
    mask_multiFour = (sum(mask_multiFour,2) == 4);
    SixCN_multiFour(i) = sum(mask_multiFour);
    mask_multiFive = (SixCNs(i,1)==FCNs |SixCNs(i,2)==FCNs |SixCNs(i,3)==FCNs |SixCNs(i,4)==FCNs |SixCNs(i,5)==FCNs |SixCNs(i,6)==FCNs);
    mask_multiFive = (sum(mask_multiFive,2) == 5);
    SixCN_multiFive(i) = sum(mask_multiFive);  
end



%% ##### 3.2.  #Corners for Ni data #####
% """
% First, get faces_An4 and faces_An5 from the first section in /Grain Tracking/main
% """
% ---- Now use all general Quads -----
QNs_An4 = result_An4{1,1};
QNs_An5 = result_An5{1,1};

C_An4 = zeros(length(faces_An4), 1);
for i = 1:length(faces_An4)
    A = faces_An4(i,1);
    B = faces_An4(i,2);
    mask_QN_A = (QNs_An4 == A);
    mask_QN_B = (QNs_An4 == B);
    mask_QN_AB = ((sum(mask_QN_A, 2) + sum(mask_QN_B, 2)) == 2);
    C_An4(i) = sum(mask_QN_AB);
end
C_An5 = zeros(length(faces_An5), 1);
for i = 1:length(faces_An5)
    A = faces_An5(i,1);
    B = faces_An5(i,2);
    mask_QN_A = (QNs_An5 == A);
    mask_QN_B = (QNs_An5 == B);
    mask_QN_AB = ((sum(mask_QN_A, 2) + sum(mask_QN_B, 2)) == 2);
    C_An5(i) = sum(mask_QN_AB);
end








