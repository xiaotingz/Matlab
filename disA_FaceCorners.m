% load data
clear

triangles_raw = textread('/Users/xiaotingzhong/Desktop/091616_STO_1470/Modified/subset1/subset1_misA_FindNagain_triangles.txt');
nodes = textread('/Users/xiaotingzhong/Desktop/091616_STO_1470/Modified/subset1/subset1_misA_FindNagain_nodes.txt');
file = ('/Users/xiaotingzhong/Desktop/091616_STO_1470/Modified/subset1/subset1_misA_FindNagain_forParaview.dream3d');
facelabel_raw = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
num_of_neigh = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
neighborList = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NeighborList'));
DisAngle = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/MisorientationList'));
surf_grain = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/SurfaceFields'));

% notice the ID of the triangles and nodes start from 0, change to start from 1
triangles_raw = triangles_raw(2:length(triangles_raw),1:4) + 1;
triangles_raw = triangles_raw.';
nodes = nodes(2:length(nodes),1:2);
nodes(:,1) = nodes(:,1) + 1;

% clear the surface triangles
tmp1 = facelabel_raw(1,:);
tmp2 = facelabel_raw(2,:);
tmp4 = triangles_raw(2,:);
tmp5 = triangles_raw(3,:);
tmp6 = triangles_raw(4,:);
boolean1 = facelabel_raw(1,:) > 0 & facelabel_raw(2,:) > 0 ;
tri_info(1,:) = tmp1(boolean1);
tri_info(2,:) = tmp2(boolean1);
tri_info(3,:) = zeros(1,length(tri_info));
tri_info(4,:) = tmp4(boolean1);
tri_info(5,:) = tmp5(boolean1);
tri_info(6,:) = tmp6(boolean1);
clear tmp1 tmp2 tmp4 tmp5 tmp6 boolean1

% prepare for count #corners
    % most faces are composed of triangles labeled [A,B] and [B,A]. Unify
    % the way of label to make faces labeled uniquely as [A,B], A<B.
tri_info = tri_info.';
for i = 1:length(tri_info)
    if tri_info(i,1) > tri_info(i,2)
        tmp = tri_info(i,1);
        tri_info(i,1) = tri_info(i,2);
        tri_info(i,2) = tmp;
    end
end
tri_info_sorted = sortrows(tri_info,[1,2]);
%%
% check every triangle on one face to count #(quad points), namely #(corners)
cnt_corner = 0;
cnt_face = 1;
labelList_use(:,1:2) = unique(tri_info(:,1:2),'rows');
face_corners = zeros(length(labelList_use),3);
for i = 1:length(tri_info_sorted)-1
%     if facelabel(i) = facelabel(i+1), the tri(i) and tri(i+1) belongs to the same face
    if tri_info_sorted(i,1) == tri_info_sorted(i+1,1) && tri_info_sorted(i,2) == tri_info_sorted(i+1,2)
        for j = 4:6
            if nodes(tri_info_sorted(i,j),2) == 4;
                cnt_corner = cnt_corner + 1;
                cornerlist(cnt_corner) = tri_info_sorted(i,j);
            end
        end
%     if not equal, tri(i) still belongs to the current face. tri(i+1) will belong to the next face
    else
        for j = 4:6
            if nodes(tri_info_sorted(i,j),2) == 4;
                cnt_corner = cnt_corner + 1;
                cornerlist(cnt_corner) = tri_info_sorted(i,j);
            end
        end
        face_corners(cnt_face,1) = tri_info_sorted(i,1);
        face_corners(cnt_face,2) = tri_info_sorted(i,2);
        face_corners(cnt_face,3) = length(unique(cornerlist));
            if face_corners(cnt_face,3) == 0
                info = ['face between [',num2str(face_corners(cnt_face,1)),', ',num2str(face_corners(cnt_face,2)), '] has zero #corners'];
                disp(info);
            end
            
        cnt_face = cnt_face + 1;
        cornerlist = [];
        cnt_corner = 0;
    end
end
%     because it take 2 triangles to do the compare, the last triangle wasn't looped
for j = 4:6
    if nodes(tri_info_sorted(i+1,j),2) == 4;
        cnt_corner = cnt_corner + 1;
        cornerlist(cnt_corner) = tri_info_sorted(i+1,j);
    end
end
face_corners(cnt_face,1) = tri_info_sorted(i+1,1);
face_corners(cnt_face,2) = tri_info_sorted(i+1,2);
face_corners(cnt_face,3) = length(unique(cornerlist));
    if face_corners(cnt_face,3) == 0
        info = ['face between [',num2str(face_corners(cnt_face,1)),', ',num2str(face_corners(cnt_face,2)), '] has zero #corners'];
        disp(info);
    end

% sort disorientation angle to match the order of face_corners
A = zeros(size(neighborList));
for i = 2:length(num_of_neigh) % 1st element in num_of_neigh is 0. Have to do this because need use i-1
    j = 1;
    while j <= num_of_neigh(i)
        cnt = sum(num_of_neigh(1:i-1)) + j;
        A(cnt) = i;
        j = j + 1;
    end
end
A = A - 1;   % grainID is matlabID-1

%%
face_disA_partial = [A, neighborList, DisAngle];
for i = 1:length(face_disA_partial)
    if face_disA_partial(i,1) > face_disA_partial(i,2)
        tmp = face_disA_partial(i,1);
        face_disA_partial(i,1) = face_disA_partial(i,2);
        face_disA_partial(i,2) = tmp;
    end
end
face_disA_partial = sortrows(face_disA_partial, [1,2]);

face_disA = zeros(length(face_disA_partial)/2,3);
cnt_uf = 0;
for i = 1:length(face_disA_partial)-1
    if face_disA_partial(i,1)==face_disA_partial(i+1,1) && face_disA_partial(i,2)==face_disA_partial(i+1,2)
        cnt_uf = cnt_uf + 1;
        face_disA(cnt_uf,1) = face_disA_partial(i,1);
        face_disA(cnt_uf,2) = face_disA_partial(i,2);
        face_disA(cnt_uf,3) = face_disA_partial(i,3);
            if roundn(face_disA_partial(i,3),-2) ~= roundn(face_disA_partial(i+1,3),-2)
                info = [num2str(i), ' th two partial face disAs not equal'];
                disp(info)
            end
    else
%         for disA list, the case that a face is labeled only once doesn't exist
        continue
    end
end


%%
figure(1)
scatter(face_disA(:,3), face_corners(:,3),'.');
set(gca,'fontsize',19)
xlabel('Disorientation Angle, °','FontSize',21);
ylabel('#Corners per Face','FontSize',21);
xlim([0,70]);

figure(2)
data_aveplot = zeros(70,3);
data_aveplot(:,1) = (1:70);
for i = 1:length(face_disA)
    data_aveplot(round(face_disA(i,3)),2) = data_aveplot(round(face_disA(i,3)),2) + face_corners(i,3);
    data_aveplot(round(face_disA(i,3)),3) = data_aveplot(round(face_disA(i,3)),3) + 1;
end
scatter(data_aveplot(:,1),data_aveplot(:,2)./data_aveplot(:,3),'filled');
set(gca,'fontsize',19)
xlabel('Disorientation Angle, °','FontSize',21);
ylabel('Average Corners per Face','FontSize',21);
line([0,70],[5.14,5.14],'LineStyle','--','color','k');
xlim([0,70]);


%% if a face was between two surface grains, it wasn't complete, delete
bool_inbulk = ones(length(face_disA),1);
for i = 1:length(face_disA)
%     the indexes in surface_grains are grainID+1
    gA = face_disA(i,1) + 1;
    gB = face_disA(i,2) + 1;
    if surf_grain(gA) == 1 && surf_grain(gB) == 1
        bool_inbulk(i) = 0;
    end
end
bool_inbulk = logical(bool_inbulk);

tmp = face_disA(:,3);
tmp2 = face_corners(:,3);
data_plot(:,1) = tmp(bool_inbulk);
data_plot(:,2) = tmp2(bool_inbulk);
clear tmp tmp2 
%%
figure(3)
scatter(data_plot(:,1), data_plot(:,2),'.');
set(gca,'fontsize',19)
xlabel('Disorientation Angle, °','FontSize',21);
ylabel('#Corners per Face','FontSize',21);
xlim([0,70]);

figure(4)
data_aveplot = zeros(70,3);
data_aveplot(:,1) = (1:70);
for i = 1:length(data_plot)
    data_aveplot(round(data_plot(i,1)),2) = data_aveplot(round(data_plot(i,1)),2) + data_plot(i,2);
    data_aveplot(round(data_plot(i,1)),3) = data_aveplot(round(data_plot(i,1)),3) + 1;
end
scatter(data_aveplot(:,1),data_aveplot(:,2)./data_aveplot(:,3),'filled');
set(gca,'fontsize',19)
xlabel('Disorientation Angle, °','FontSize',21);
ylabel('Average Corners per Face','FontSize',21);
line([0,70],[5.14,5.14],'LineStyle','--','color','k');
xlim([0,70]);
       
        
        
        
        


