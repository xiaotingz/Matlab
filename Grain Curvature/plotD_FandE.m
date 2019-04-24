
file = '/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2_GBCD.dream3d';
% file = ('/Users/xiaotingzhong/Desktop/Datas/Jan.31 Austenite/180402_austenite_recons.dream3d');
% --------------------------- V6 structure ---------------------------
numNeigh = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
NeighborList = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
grain_diameter_raw = double(roundn(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters'),-5)).';
surfGrain = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')).';
numNeigh(1) = [];
grain_diameter_raw(1) = [];
surfGrain(1) = [];
numEdges = calcNumEdges(numNeigh, NeighborList);



% load('STOcombinedCentroThres_grain_ForCal.mat');
% grain_ForCal = [grainId, grainDiameter, #Faces, #edges]
sizes = grain_diameter_raw(logical(surfGrain));
numFaces = numNeigh(logical(surfGrain));
numEdges = numEdges(logical(surfGrain));


figure(1)
scatter(sizes, numEdges/sum(numEdges)*length(numEdges),'filled','r')
hold on
scatter(sizes, numFaces/sum(numFaces)*length(numFaces),'filled','b')
xlabel('D, \mum','FontSize',21);
ylabel('Normalized Topology Parameter','FontSize',21);
legend({'# Faces / <# Faces>','# Edges / <# Edges>'},'FontSize',17,'Location','northwest');
 
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

figure(2)
scatter(numEdges, numFaces, 'filled', 'k')
xlabel('# Faces','FontSize',21);
ylabel('# Edges','FontSize',21);

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