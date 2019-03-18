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
%             output: grain_ForCal = [grainId, grainDiameter, #Faces, #edges]
%         calcGrainCurvature
%             objective: assemble the integral face curvatures for integral grain curvature
%             output: data_grain = [grainId, grainDiameter, #Faces, #edges, IntegralGrainCurvature];
%         gridData
%             objective: bin the grain data and plot
%             output: data_grid = [binStart, binEnd, s_x/cnt, s_cur/cnt, cnt]
% - If multiple volumes
%         1. calc data_grain for the different volumes
%         2. data_grain = [data_grain_sub1; data_grain_sub2]; 
%         3. plot with the concatenated data_grain
% ##################################################################
clear

% subset1
% file = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2_GBCD.dream3d');
% centro_file = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2_GBCD.dream3d');
file = ('/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/problemData/Jan31_Aa0_CurvDistri10.dream3d');
centro_file = ('/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/problemData/Jan31_Aa0_CurvDistri10.dream3d');
X=434*0.15; Y=267*0.15; Z=100*0.2;
% subset2
% file = ('/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/Jan31_Fca0.dream3d');
% centro_file = ('/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/Jan31_Fca0.dream3d');
% X=234*0.15; Y=267*0.15; Z=68*0.2;
% -- criterions in filterGrains: 'centroidPos' | 'touchingFS' | 'numFaces'| 'NN_centoridPos' | 'NN_touchingFS'
criterion = 'NN_centoridPos';
% -- xAxis in gridData: 'D' | 'numFaces' | 'numEdges'
xAxis = 'numFaces';

%  -------------------------- bounding boxes -------------------------- 
% Austenite --- X=434*0.15; Y=267*0.15; Z=100*0.2;
% Ferrite   --- X=234*0.15; Y=267*0.15; Z=68*0.2;
% STO_1470_sub1 --- X=232*0.3; Y=129*0.3; Z=36*0.3;
% STO_1470_sub2 --- X=213*0.3; Y=297*0.3; Z=40*0.3;
% Mg_undeform --- X=167*0.15; Y=133*0.15; Z=73*0.15;
% Mg_small --- X=167*0.15; Y=133*0.15; Z=17*0.15;
% MNK_Ti --- X=431*0.5; Y=109*0.5; Z=104*0.3;

% -------------------------- load data_v4 -------------------------- 
centroids = roundn(h5read(centro_file,'/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';
grain_diameter_raw = roundn(h5read(centro_file,'/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
facelabel = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
curvature_of_triangle = h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures').';
num_of_neigh = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
neighborList = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NeighborList'));
triangle_area_raw = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-8).';
%  -------------------------- load data_v6 -------------------------- 
% centroids = abs(roundn(h5read(centro_file,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).');
% grain_diameter_raw = roundn(h5read(centro_file,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters'),-5).';
% num_of_neigh = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
% neighborList = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
% facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'));
% curvature_of_triangle = h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures');
% triangle_area_raw = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-8);

centroids(1,:) = [];
grain_diameter_raw(1) = [];
num_of_neigh(1) = [];
data_raw = [facelabel; curvature_of_triangle; triangle_area_raw];


%  -------------------------- get the grid bin data -------------------------- 
data_face = calcFaceCurvature(data_raw);
grain_ForCal = filterGrains(criterion, facelabel, num_of_neigh, neighborList, X,Y,Z, centroids, grain_diameter_raw);
data_grain = calcGrainCurvature(data_face, grain_ForCal);


%% -------------------------- make G from Ms -------------------------- 
% % data_grain = [grainId, grainDiameter, #Faces, #edges, IntegralGrainCurvature];
% data_grain(:,5) = -data_grain(:,5)./(data_grain(:,2).*((4/3*pi)^(1/3)/2));
%
%  -------------------------- plot average bin data & the standard deviations -------------------------- 
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

%  -------------------------- plot parameters -------------------------- 
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
elseif strcmp (xAxis, 'numEdges')
    ax = gca;
    % ax.XTick = [0:5:40];
    xlabel('E','FontSize',21);
    xlim([0,100]);
    ylim([-30,10]);
    % text(0.55,-26,'(a)','FontWeight','bold','FontSize',30)
    line([0,max(data_grain(:,4))+50],[0,0],'LineStyle','--', 'Color',[0.5 0.5 0.5])
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
% % -------------------------------------------------------------------------------

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