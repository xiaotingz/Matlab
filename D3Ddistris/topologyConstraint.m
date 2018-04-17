% ############################################################################
% Topology constraint for the entire sample
% NOTES
% - for 3D structure, - C + F - E + V = 1
% ############################################################################


file = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2_GBCD.dream3d');
% file = ('/Users/xiaotingzhong/Desktop/Datas/Jan.31 Austenite/180402_austenite_recons.dream3d');
% triNodes = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList')).';
nodeTypes = double(h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType'));
facelabel_raw = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
numNeigh = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors'));
NeighborList = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
surfGrain = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures'));
numNeigh(1) = [];
surfGrain(1) = [];
grainIDs = [1:1:length(numNeigh)];
surfGrainID = grainIDs(logical(surfGrain));
innerGrainID = grainIDs(~logical(surfGrain));

C = length(numNeigh) - 1;

E_inside = sum(numEdges)/3;  
edgesFreeSurf = [];
for i = 1:length(surfGrainID)
    list1 = getNeighList(surfGrainID(i), numNeigh, NeighborList);
    intersection_inner = intersect(list1, surfGrainID);
    edgesFreeSurf = vertcat(edgesFreeSurf, [ones(length(intersection_inner),1)*surfGrainID(i), intersection_inner]);
end
tmp = sort(edgesFreeSurf,2);
unique(tmp,'rows');
if length(unique(tmp,'rows')) ~= length(edgesFreeSurf)/2
    warning('Hey the number of edges on free surface is NOT RIGHT!');
end
E_freeSurf = length(edgesFreeSurf)/2;
E_all = E_inside + E_freeSurf;


F_inside = sum(numNeigh)/2;
mask = (all(facelabel_raw > 0, 2));
surfTris = facelabel_raw(~mask,:);
F_surface = length(unique(surfTris, 'rows'));
F_all = F_inside + F_surface;

mask = (nodeTypes == 4);
V_inside = sum(mask);
mask = (nodeTypes == 14);
V_surface = sum(mask);
V_all = V_inside + V_surface;

- C + F_inside - E_inside + V_inside


- C + F_all - E_all + V_all

%%  the statistics for the thresholded grains

numEdges_thres = numEdges(~logical(surfGrain.'));
numNeigh_thres = numNeigh(~logical(surfGrain)).';

% -- #Grains & #Vertices are simple
C_thres = sum(1-surfGrain);
V_thres = V_inside;

% -- #Faces = #InnerFaces/2 + #BorderFaces
numBorderFaces = 0;
% numInnerFaces = 0;
for i = 1:length(innerGrainID)
    numBorderFaces = numBorderFaces + length(intersect(getNeighList(innerGrainID(i), numNeigh, NeighborList), surfGrainID));
%     numInnerFaces = numInnerFaces + length(intersect(getNeighList(innerGrainID(i), numNeigh, NeighborList), innerGrainID));
end
F_thres = (sum(numNeigh_thres) - numBorderFaces)/2 + numBorderFaces;


% -- #Edges = #InnerEdges/3 + #BetweenInnerAndBorderEdges/2 + #BorderEdges
edgeShare3 = [];
edgeShare2 = [];
edgeShare1 = [];
for i = 1:length(numNeigh_thres)
% for i = 1458:1458
    keyGrain = innerGrainID(i);
    keyNeighbors = getNeighList(keyGrain, numNeigh, NeighborList);

%     commonNeighPairs = [];
    for j = 1:numNeigh(keyGrain)
        toLook = keyNeighbors(j);
        list2 = getNeighList(toLook, numNeigh, NeighborList);
        intersection = Intersect(keyNeighbors, list2);
        for k = 1:length(intersection)
            if ismember(toLook, innerGrainID) && ismember(intersection(k), innerGrainID)
                edgeShare3 = vertcat(edgeShare3, intersection(k));
            elseif isempty(intersect([toLook,intersection(k)], innerGrainID))
                edgeShare1 = vertcat(edgeShare1, intersection(k));
            else
                edgeShare2 = vertcat(edgeShare2, intersection(k));
            end  
        end
    end

end
if (length(edgeShare1) + length(edgeShare2) + length(edgeShare3))/2 ~= sum(numEdges_thres)
    warning('Hey the #Edges for the inner grains are NOT RIGHT!');
end


E_thres = length(edgeShare3)/3 + length(edgeShare2)/2 + length(edgeShare1);
% 
- C_thres + F_thres - E_thres + V_thres
% 
% 


%% 
% ############################################################################
% Topology constraint for a single grain
% NOTES
% - for an inside grain, F - E + V = 2
% - for a grain touching the free surface, F - E + V = 1
% ############################################################################
file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_mesh.dream3d');
% file = ('/Users/xiaotingzhong/Desktop/Datas/Jan.31 Austenite/180402_austenite_recons.dream3d');
nodeTypes = double(h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType'));
TriNodes = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')).';
facelabel_raw = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
numNeigh = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors2'));
NeighborList = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList2'));
surfGrain = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures'));
numNeigh(1) = [];
surfGrain(1) = [];
grainIDs = [1:1:length(numNeigh)];
surfGrainID = grainIDs(logical(surfGrain));
innerGrainID = grainIDs(~logical(surfGrain));
TriNodes = TriNodes + 1;

mask = (all(facelabel_raw > 0, 2));
innerTris = facelabel_raw(mask, :);
innerTriNodes = TriNodes(mask, :);

TopologyInfo =[];
for i = 1:length(innerGrainID)
    objGrain = innerGrainID(i);
    mask = (innerTris(:,1) == objGrain | innerTris(:,2) == objGrain);
%     objGrainTris = innerTris(mask,:);
    objGrainTrisNodes = innerTriNodes(mask,:);
    mask2 = (nodeTypes(objGrainTrisNodes) == 4);
    V = length(unique(objGrainTrisNodes(mask2)));
    F = numNeigh(objGrain);
    E = numEdges(objGrain);
    TopologyInfo = vertcat(TopologyInfo, [objGrain, V, E, F]);
end

tmp = numNeigh(~logical(surfGrain));
getNeighList(objGrain, numNeigh, NeighborList)

NNinfo = [grainIDs(TopologyInfo(:,1)).', numNeigh(TopologyInfo(:,1)).', TopologyInfo(:,2) - TopologyInfo(:,3) + TopologyInfo(:,4)];


histogram(NNinfo(:,3));
title('Topology of Inner Grains, Ni\_An4','FontSize',21)
xlabel('V - E + F','FontSize',21);
ylabel('# Grains','FontSize',21);
ax = gca;
set(ax,'fontsize',19)
set(gca,'FontWeight','bold','linewidth',2)
% disable tick on top and right of the box
    % get handle to current axes
a = gca;
    % set box property to off and remove background color
set(a,'box','off','color','none')
    % create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'FontWeight','bold','linewidth',2);
    % set original axes as active
axes(a)
    % link axes in case of zooming
linkaxes([a b])





% tmp2 = (facelabel_raw(:,1) == objGrain | facelabel_raw(:,2) == objGrain);
% tmp2 = facelabel_raw(tmp2, :);
% unique(tmp2, 'rows')







