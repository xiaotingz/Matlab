% load data
file = ('/Volumes/RESEARCH/Grain Curvature/MinSize=100/setTo0/Jan31_A_setTo0.dream3d');

num_of_neigh = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
neighborList = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NeighborList'));
DisAngle = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/MisorientationList'));
% triangle_area_raw = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-8);
facelabel = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
triCurv = abs(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));
data_raw = [facelabel; triCurv.'];

tmp1 = data_raw(1,:);
tmp2 = data_raw(2,:);
tmp3 = data_raw(3,:);

% get rid of bad datas: 
    % 1.facelabel <= 0; 2.extreme curvature value; 3.NaN(sometimes)
boolean1 = data_raw(1,:) > 0 & data_raw(2,:) > 0 & ~isnan(data_raw(3,:)) & data_raw(3,:) < 100 & data_raw(3,:) > -100;
data_cleared(1,:) = tmp1(boolean1);
data_cleared(2,:) = tmp2(boolean1);
data_cleared(3,:) = tmp3(boolean1);

data_sorted = sortrows(data_cleared.');

% create lableList[A,B] from num_of_neigh and neighborList(B)
A = zeros(size(neighborList));
for i = 2:length(num_of_neigh) % fisrt element in num_of_neigh is 0
    j = 1;
    while j <= num_of_neigh(i)
        cnt = sum(num_of_neigh(1:i-1)) + j;
        A(cnt) = i;
        j = j + 1;
    end
end
A = A - 1;   % grainID is matlabID-1

% labelList in use is the unique labels of triangles
labelList_full = [A, neighborList];
labelList_use = unique(data_sorted(:,1:2),'rows');
disA_use = zeros(length(labelList_use),1);

cnt = 1;
for i = 1:length(labelList_full)
    if labelList_full(i,1) == labelList_use(cnt,1) && labelList_full(i,2) == labelList_use(cnt,2) 
        disA_use(cnt) = DisAngle(i);
        cnt = cnt + 1;
    else
        continue
    end
end
clear data_raw boolean1 data_cleared file tmp1 tmp2 tmp3 A

%% 
% for every face in labelList_use, sum curvature, assign misorientation
    % data_final = [#triangles on the face, total curvature on the face, disorientation]
data_final = zeros(length(labelList_use),3);
j = 1;  % j -- jth faces
k = 1;   % k --- kth traingle on the jth face
total = 0;   % total --- total curvature
for i = 1:length(data_sorted)-1
        % compare if label of kth triangle equals (k+1)th, if yes, sum curvature
    if data_sorted(i,1) == data_sorted(i+1,1) && data_sorted(i,2) == data_sorted(i+1,2)
        total = total + data_sorted(i,3);
        k = k + 1;
        % if not, sum kth into this face, start a new face
    else 
        total = total + data_sorted(i,3);
        data_final(j,1) = k;
        data_final(j,2) = total;
        data_final(j,3) = disA_use(j);
        j = j + 1;
        k = 1;
        total = 0;
    end
end
data_final(j,2) = total + data_sorted(i+1,3);
data_final(j,1) = k;
data_final(j,3) = disA_use(j);

%%
% sort curvatures according to misorientation
start = 0;
nstep = 35;
stepsize = 2;
data_grid = zeros(nstep,2);
totalCurv = 0;
totalNtri = 0;
for j = 1 : nstep
    for i = 1 : length(data_final)
        % data_final = [#triangles on the face, total curvature on the face, disorientation]
        if data_final(i,3) > start && data_final(i,3) <= start + stepsize
            totalCurv = totalCurv + data_final(i,2);
            totalNtri = totalNtri + data_final(i,1);
        end
    end
    data_grid(j,1) = start + stepsize/2;
    data_grid(j,2) = totalCurv/totalNtri;
    start = start + stepsize;
    totalCurv = 0;
    totalNtri = 0;
end

%% 
figure(1)
h1 = histogram(DisAngle);
h1.Normalization = 'probability';
h1.BinWidth = 2;
set(h1,'facecolor',[0.5,0.5,0.5]);
set(gca,'fontsize',14)
xlabel('Disorientation Angle, °','FontSize',17);
ylabel('Frequency','FontSize',17);
% disable tick on top and right of the box
a = gca;set(a,'box','off','color','none');
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
axes(a)
linkaxes([a b])


figure(2)
scatter(data_grid(:,1), data_grid(:,2),'filled','k'); hold on
% scatter(60,0.3,90,'filled','p');
% text(36,0.32,'coherent twin boundary \rightarrow','FontSize',12);
box on
set(gca,'fontsize',14)
xlabel('Disorientation Angle, °','FontSize',17);
ylabel('Average Curvature, \mum^{-1}','FontSize',17);
% disable tick on top and right of the box
a = gca;set(a,'box','off','color','none');
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
axes(a)
linkaxes([a b])