% ##################################################################
% NOTES:
%  - Variable to change:
%         File name
%         Bounding box: X, Y, Z
%         criterion: the criterion chosing the grains to calculate
%         xAxis: the xAxis for plot 
%  - Functions:
%         calcFaceCurvature
%             objective: calculate the integral curvature for each face
%             output: data_face[labelA, labelB, faceArea, faceIntegralCurv/faceArea]
%         filterGrains
%             objective: exclude the free surface grains
%             output: grain_ForCal = [grainId, grainDiameter, #Faces]
%         calcGrainCurvature
%             objective: assemble the integral face curvatures for integral grain curvature
%             output: data_grain = [grainId, grainDiameter, #Faces, IntegralGrainCurvature];
%         gridData
%             objective: bin the grain data and plot
%             output: data_grid = [binStart, binEnd, s_x/cnt, s_cur/cnt, cnt]
% ##################################################################
clear
file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_mesh.dream3d');
% --- xAxis is used for plot, can be either 'numFaces'|'D' ---
xAxis = 'D';

%  -------------------------- load data -------------------------- 
centroids = abs(roundn(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).');
grain_diameter = roundn(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters'),-5).';
num_of_neigh = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
neighborList = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'));
tri_curv = h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures');
tri_area = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-8);

centroids(1,:) = [];
grain_diameter(1) = [];
num_of_neigh(1) = [];
data_raw = [facelabel; tri_curv; tri_area];


%  -------------------------- grain curvature data -------------------------- 
data_face = calcFaceCurvature(data_raw);
grain_id = (1:length(num_of_neigh))';
% """
% - One can filter the grains according to different criterions. However
% that filter script was written for cubic shape samples so you probably
% want to modify it before using. 
% grain_ForCal = filterGrains(criterion, facelabel, num_of_neigh, neighborList, X,Y,Z, centroids, grain_diameter_raw);
% - Note there were some problem in the calculation for number of edges. I
% haven't corrected that. The data is put here just for format needs. 
% """
numEdges = calcNumEdges(num_of_neigh, neighborList);
grain_ForCal = [grain_id, grain_diameter, num_of_neigh, numEdges];
data_grain = calcGrainCurvature(data_face, grain_ForCal);


% %% -------------------------- make G from Ms -------------------------- 
% data_grain(:,5) = -data_grain(:,5)./(data_grain(:,2).*((4/3*pi)^(1/3)/2));


%  -------------------------- plot average bin data & the standard deviations -------------------------- 
% """
% To change the grid parameters of grids, go into gridData and change the
% width paramter.
% """
data_grid = gridData(xAxis, data_grain);
scatter(data_grid(:,3),data_grid(:,4),50,'filled','o','k');
% -- plot line: standard mean deviation
range = zeros(length(data_grid),2);
for i = 1:length(data_grid)
    range(i,1) = data_grid(i,4) + data_grid(i,6); % up value
    range(i,2) = data_grid(i,4) - data_grid(i,6); % down value
    line([data_grid(i,3),data_grid(i,3)],[range(i,1),range(i,2)],'color','k');
    hold on
end

%  -------------------------- plot axis -------------------------- 
if strcmp (xAxis, 'numFaces')
    ax = gca;
    % ax.XTick = [0:5:40];
    xlabel('F','FontSize',21);
%     xlim([0,40]);
%     ylim([-40,10]);
    % text(0.55,-26,'(a)','FontWeight','bold','FontSize',30)
    line([0,max(data_grain(:,3))+50],[0,0],'LineStyle','--', 'Color',[0.5 0.5 0.5])
elseif strcmp (xAxis, 'D')
    ax = gca;
    % ax.XTick = [0:5:40];
    xlabel('D (\mum)','FontSize',21);
    % text(0.55,-26,'(a)','FontWeight','bold','FontSize',30)
    line([0,max(data_grain(:,2))+3],[0,0],'LineStyle','--', 'Color',[0.5 0.5 0.5])
% elseif strcmp (xAxis, 'numEdges')
%     ax = gca;
%     % ax.XTick = [0:5:40];
%     xlabel('E','FontSize',21);
%     xlim([0,100]);
%     ylim([-30,10]);
%     % text(0.55,-26,'(a)','FontWeight','bold','FontSize',30)
%     line([0,max(data_grain(:,4))+50],[0,0],'LineStyle','--', 'Color',[0.5 0.5 0.5])
end

% % ---------------------------- G from Ms, fit linear line -----------------------
% xlim([0,40])
% ylabel('$\mathcal{G''}$','Interpreter','latex','FontSize',21);
% x = data_grid(:,3);
% y = data_grid(:,4);
% mask = (~isnan(x) & ~isnan(y)); 
% x = x(mask);
% y = y(mask);
% X = [ones(length(x),1) x];
% b = X\y;
% yCalc2 = X*b;
% plot(x,yCalc2,'-.','LineStyle','-.', 'Color',[0.5 0.5 0.5])

% ---------------------------- adjust plot details -----------------------
set(ax,'fontsize',19)
ylabel('M_{S} (\mum)','FontSize',21);
% set(gca,'FontWeight','bold','linewidth',2)
% disable tick on top and right of the box
    % get handle to current axes
a = gca;
    % set box property to off and remove background color
set(a,'box','off','color','none')
    % create new, empty axes with box but without ticks
% b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'FontWeight','bold','linewidth',2);
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
% set original axes as active
axes(a)
    % link axes in case of zooming
linkaxes([a b])