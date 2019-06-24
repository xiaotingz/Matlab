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
%             output: data_grid[binStart, binEnd, s_x/cnt, s_cur/cnt, cnt]
% - If multiple volumes
%         1. calc data_grain for the different volumes
%         2. data_grain = [data_grain_sub1; data_grain_sub2]; 
%         3. plot with the concatenated data_grain
% ##################################################################
% clear

% % file = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub1_GBCD.dream3d');
% % centro_file = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub1_GBCD.dream3d');
X=398*2.8125; Y=404*2.8125; Z=104*4.0;
% -- criterions in filterGrains: 'centroidPos' | 'touchingFS' |
% 'numFaces' | 'NN_centoridPos' | 'NN_touchingFS'
criterion = 'none';
% -- data_grid limits
start = -80;
width = 2;
step = 80;

%  -------------------------- bounding boxes -------------------------- 
% Austenite --- X=434*0.15; Y=267*0.15; Z=100*0.2;
% Ferrite   --- X=234*0.15; Y=267*0.15; Z=68*0.2;
% STO_1470_sub1 --- X=232*0.3; Y=129*0.3; Z=36*0.3;
% STO_1470_sub2 --- X=213*0.3; Y=297*0.3; Z=40*0.3;
% Mg_undeform --- X=167*0.15; Y=133*0.15; Z=73*0.15;
% Mg_small --- X=167*0.15; Y=133*0.15; Z=17*0.15;
% MNK_Ti --- X=431*0.5; Y=109*0.5; Z=104*0.3;

% -------------------------- load data_v4 -------------------------- 
% centroids = roundn(h5read(centro_file,'/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';
% grain_diameter_raw = roundn(h5read(centro_file,'/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
% facelabel = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
% curvature_of_triangle = h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures').';
% num_of_neigh = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
% neighborList = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NeighborList'));
% triangle_area_raw = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-8).';

file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
file_an5 = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth.dream3d';
centro_file = file_an4;
file = centro_file;
% file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
% centro_file = file_an4;
% file = file_an4;

%  -------------------------- load data_v6 -------------------------- 
centroids = abs(roundn(h5read(centro_file,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).');
grain_diameter_raw = roundn(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters'),-5).';
num_of_neigh = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
neighborList = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'));
curvature_of_triangle = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5);
triangle_area_raw = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5);
triangle_min_ang = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5);
% surfGrains = logical(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures'));

% surfGrains(1) = [];
centroids(1,:) = [];
grain_diameter_raw(1) = [];
num_of_neigh(1) = [];
data_raw = [facelabel; curvature_of_triangle; triangle_area_raw; triangle_min_ang];


%  -------------------------- get the grid bin data -------------------------- 
% data_face = calcFaceCurvature(data_raw);
eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;
data_face = calcFaceCurvature(data_raw, eps_curv, eps_area, eps_min_ang);

% --- grain_ForCal = [ID_ForCal, D_ForCal, NNeigh_ForCal, numEdges_ForCal] ---;
grain_ForCal = filterGrains(criterion, facelabel, num_of_neigh, neighborList, X,Y,Z, centroids, grain_diameter_raw);
% --- data_grain = [ID_ForCal, D_ForCal, NNeigh_ForCal, numEdges_ForCal, grain_itg_curv] ---;
data_grain = calcGrainCurvature(data_face, grain_ForCal);


ID_ForCal = data_grain(:,1);
F_mF_diff = zeros(size(ID_ForCal));
for i = 1 : length(ID_ForCal)
    grainID = ID_ForCal(i);
    total = 0;
    nlist_start = sum(num_of_neigh(1:(grainID-1))) + 1;
    nlist_end = sum(num_of_neigh(1:grainID));
        for j = nlist_start : nlist_end
            k = neighborList(j);
            total = total + num_of_neigh(k);
        end
    F_mF_diff(i) = num_of_neigh(grainID) - total/num_of_neigh(grainID);
end

curv_FmF = [F_mF_diff, data_grain(:,5)./(data_grain(:,2).*((4/3*pi)^(1/3)/2))];

% %% ------------------------------ filter STO data ------------------------------
% % load(G_FmF_STO1470_allGrains.mat)
% % load(data_grain_STO1470combined_thresCentroid.mat)
% curv_FmF_1 = curv_FmF(1:720, :);
% curv_FmF_2 = curv_FmF(721:end, :);
% curv_FmF_1 = curv_FmF_1(data_grain(1:134, 1), :);
% curv_FmF_2 = curv_FmF_2(data_grain(135:end, 1), :);
% curv_FmF = [curv_FmF_1; curv_FmF_2];

%%
scatter(curv_FmF(:,1),curv_FmF(:,2),10,[0.5,0.5,0.5],'filled');
hold on

start = -80;
width = 1;
step = 160;

data_grid = zeros(step,3);
cnt = 0;
total = 0;
for i = 1:step
% mean Curv for each Bin
    for j = 1:length(curv_FmF)
        if curv_FmF(j,1) >= start && curv_FmF(j,1) < start+width
            total = total + curv_FmF(j,2);
            cnt = cnt + 1;
        end
    end
    data_grid(i,1) = start + width/2;
    data_grid(i,2) = total/cnt;
    data_grid(i,3) = cnt;
    
    start = start + width;
    total = 0;
    cnt = 0;
end
data_grid(data_grid(:,3) <= 3,:) = [];
scatter(data_grid(:,1), data_grid(:,2), 60, 'r', 'filled', 's')
% start = -20;
% width = 1;
% step = 40;
% xlim([-50, 50])
% ylim([-10, 10])

xlabel('F - <F_{NN}>','FontSize',21)
ylabel('$\mathcal{G''}$','Interpreter','latex','FontSize',21);
% ylabel('$\Delta V$','Interpreter','latex','FontSize',21);
line([-60,60],[0,0],'LineStyle','--', 'Color',[0.5 0.5 0.5])
line([0,0],[-20,20],'LineStyle','--', 'Color',[0.5 0.5 0.5])


set(gca,'fontsize',19)
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
