% test the correspondence of triangles. Find surface triangles: bool1 from
% nodes, bool2 from faceLabel, check if they are the same
clear

triangles_raw = textread('/Users/xiaotingzhong/Desktop/Datas/Jan.31 Austenite/Jan31_A_croptriangle.txt');
nodes = textread('/Users/xiaotingzhong/Desktop/Datas/Jan.31 Austenite/Jan31_A_cropnode.txt');
file = ('/Users/xiaotingzhong/Desktop/Datas/Jan.31 Austenite/Jan31_A_crop_misA.dream3d');
facelabel_raw = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));

triangles_raw = triangles_raw(2:length(triangles_raw),1:4) + 1;
nodes = nodes(2:length(nodes),1:2);
nodes(:,1) = nodes(:,1) + 1;

boolean1 = ones(length(triangles_raw),1);

for i = 1 : length(triangles_raw)
    for j = 2:4
        if nodes(triangles_raw(i,j),2) > 10
            boolean1(i) = 0;
        end
    end
end

boolean2 = facelabel_raw(1,:) > 0 & facelabel_raw(2,:) > 0 ;
boolean2 = double(boolean2.');

cnt_quad = 1;
for i = 1:length(boolean1)
    if boolean1(i) ~= boolean2(i)
        table_1minus2(cnt_quad,1) = i;
        table_1minus2(cnt_quad,2) = boolean1(i) - boolean2(i);
        cnt_quad = cnt_quad + 1;
    end
end

table_nodes = zeros(length(table_1minus2),5);
for i = 1:length(table_1minus2)
    table_nodes(i,1) = table_1minus2(i,1);
    table_nodes(i,2) = 0;
    table_nodes(i,3) = nodes(triangles_raw(table_1minus2(i,1),2),2);
    table_nodes(i,4) = nodes(triangles_raw(table_1minus2(i,1),3),2);
    table_nodes(i,5) = nodes(triangles_raw(table_1minus2(i,1),4),2);
    if (table_nodes(i,3) + table_nodes(i,4) + table_nodes(i,5)) < 10 or (table_nodes(i,3) + table_nodes(i,4) + table_nodes(i,5)) > 30
        disp('not right!')
    end
end


for i = 1:length(boolean2)
    if boolean2(i) == 0 
        surftri_nodes1 = nodes(triangles_raw(i,2),2);
        surftri_nodes2 = nodes(triangles_raw(i,3),2);
        surftri_nodes3 = nodes(triangles_raw(i,4),2);
        if surftri_nodes1 + surftri_nodes2 + surftri_nodes3 < 30
            disp('not right!')
        end
    end
end


%%
bool = (face_corners(:,1:2) == labelList_use);
sum(bool)
    
%%
clear FaceQuads FaceQuads_unique FaceQuats_coor FaceTris FaceTris_uniqu FaceTris_coor

gA = 217; gB = 688;
cnt_quad = 0;
cnt_tri = 0;
for i = 1:length(facelabel_raw)
    if (facelabel_raw(1,i)==gA && facelabel_raw(2,i)==gB) || (facelabel_raw(1,i)==gB && facelabel_raw(2,i)==gA)
        for j = 2:4
            if nodes(triangles_raw(j,i),2) == 4
                cnt_quad = cnt_quad + 1;
                FaceQuads(cnt_quad) = triangles_raw(j,i);
            end
            if nodes(triangles_raw(j,i),2) == 3
                cnt_tri = cnt_tri + 1;
                FaceTris(cnt_tri) = triangles_raw(j,i);
            end
        end
    end
end
FaceQuads_unique = unique(FaceQuads);
FaceTris_unique = unique(FaceTris);

% nodes = textread('/Users/xiaotingzhong/Desktop/091616_STO_1470/Modified/subset1/subset1_misA_FindNagain_nodes.txt');
% nodes = nodes(2:length(nodes),:);
% nodes(:,1) = nodes(:,1) + 1;

FaceQuats_coor = zeros(length(FaceQuads_unique),5);
for i = 1:length(FaceQuads_unique)
    FaceQuats_coor(i,:) = nodes(FaceQuads_unique(i),:);
end
sortrows(FaceQuats_coor,3);
scatter3(FaceQuats_coor(:,3),FaceQuats_coor(:,4),FaceQuats_coor(:,5),'filled','r')
hold on

FaceTris_coor = zeros(length(FaceTris_unique),5);
for i = 1:length(FaceTris_unique)
    FaceTris_coor(i,:) = nodes(FaceTris_unique(i),:);
end
scatter3(FaceTris_coor(:,3),FaceTris_coor(:,4),FaceTris_coor(:,5),'filled','k')




