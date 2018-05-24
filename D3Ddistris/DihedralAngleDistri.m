% 1. Read in data
% % Austenite
% tri_full = textread('/Users/xiaotingzhong/Desktop/for_FindDAD/Jan31_A_triangle.txt');
% node_full = textread('/Users/xiaotingzhong/Desktop/for_FindDAD/Jan31_A_node.txt');
% normal_full = roundn(h5read('/Users/xiaotingzhong/Desktop/for_FindDAD/Jan31_A_stats.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals'),-5);
% % Ferrite
% tri_full = textread('/Volumes/DATAS/Datas/Ferrite/Jan.31 Ferrite/Jan31_F_triangle.txt');
% node_full = textread('/Volumes/DATAS/Datas/Ferrite/Jan.31 Ferrite/Jan31_F_node.txt');
% normal_full = roundn(h5read('/Volumes/DATAS/Datas/Ferrite/Jan.31 Ferrite/Jan31_Fc_stats.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals'),-5);
% Cropped
% tri_full = textread('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_A_croptriangle.txt');
% node_full = textread('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_A_cropnode.txt');
% normal_full = roundn(h5read('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_A_cropstats.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals'),-5);
% % Nickel
tri_full = textread('/Users/xiaotingzhong/Desktop/Datas/060216 Ni/060216_Ni_triange.txt');
node_full = textread('/Users/xiaotingzhong/Desktop/Datas/060216 Ni/060216_Ni_node.txt');
normal_full = roundn(h5read('/Users/xiaotingzhong/Desktop/Datas/060216 Ni/060216_Ni_stats.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals'),-5);
% % Nickel_2
% tri_full = textread('/Volumes/DATAS/Datas/060216 Ni_2/060216_Ni_2_triangle.txt');
% node_full = textread('/Volumes/DATAS/Datas/060216 Ni_2/060216_Ni_2_node.txt');
% normal_full = roundn(h5read('/Volumes/DATAS/Datas/060216 Ni_2/060216_Ni_2_stats.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals'),-5);
% % Nickel_1
% tri_full = textread('/Users/xiaotingzhong/Desktop/060216_Ni_1/060216_Ni_1_triangle.txt');
% node_full = textread('/Users/xiaotingzhong/Desktop/060216_Ni_1/060216_Ni_1_node.txt');
% normal_full = roundn(h5read('/Users/xiaotingzhong/Desktop/060216_Ni_1/060216_Ni_1_stats.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals'),-5);


    % adjust count starting from 0 to 1
node_full(1,:) = [];
node_full(:,1) = node_full(:,1) + 1;
tri_full(1,:) = [];
tri_full(:,1:4) = tri_full(:,1:4) + 1;
tri_full(:,5:7) = [];
normal_full = normal_full.';

% 2. Find the triangles on tripleline
    % first clean surface triangles
tmp1 = tri_full(:,1); tmp2 = tri_full(:,2); tmp3 = tri_full(:,3); 
tmp4 = tri_full(:,4); tmp5 = tri_full(:,5); tmp6 = tri_full(:,6);
bool1 = (tri_full(:,5) > 0 & tri_full(:,6) > 0);
GBtri(:,1) = tmp1(bool1); GBtri(:,2) = tmp2(bool1); GBtri(:,3) = tmp3(bool1);
GBtri(:,4) = tmp4(bool1); GBtri(:,5) = tmp5(bool1); GBtri(:,6) = tmp6(bool1);
clear tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 bool1

    % then check node type, if 2 of the three nodes are tripleline node (type3), keep
bool2 = zeros(length(GBtri),1);
for i = 1:length(GBtri)
    cnt = 0;
    for j = 2:4
        nodeType = node_full(GBtri(i,j),2);
        if nodeType == 3
            cnt = cnt + 1;
        end
    end
    
    if cnt == 2
        bool2(i) = 1;
    end
end
bool2 = logical(bool2);
tmp1 = GBtri(:,1); tmp2 = GBtri(:,2); tmp3 = GBtri(:,3); 
tmp4 = GBtri(:,4); tmp5 = GBtri(:,5); tmp6 = GBtri(:,6);
TLtri(:,1) = tmp1(bool2); TLtri(:,2) = tmp2(bool2); TLtri(:,3) = tmp3(bool2);
TLtri(:,4) = tmp4(bool2); TLtri(:,5) = tmp5(bool2); TLtri(:,6) = tmp6(bool2);
clear tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 bool2 
    
% 3. Sort the triangles into groups of 3
    % creat a Tnode list to record the 2 triple line nodes of each
    % triangle, when find group, compare only these two nodes
TLnodelist = zeros(length(TLtri),2);
for i = 1:length(TLtri)
    cnt = 1;
    for j = 2:4
        nodeID = TLtri(i,j);
        nodeType = node_full(nodeID,2);
        if nodeType == 3
            TLnodelist(i,cnt) = nodeID;
            cnt = cnt + 1;
        end
    end
end
    % group triangles according to the 2 common triple line node
        % sort the two nodes to make compare convenient
        % first within a triangle(row), then through the list(colume)
TLnodelist = sort(TLnodelist,2);
TLnodelist = [TLtri(:,1),zeros(length(TLtri),1),TLnodelist];
TLnodelist = sortrows(TLnodelist,3);
    % if the group contain more or less than 3 triangles, ignore
bool_groupof3 = zeros(length(TLnodelist),1);
cnt = 1;
for i = 1:(length(TLnodelist)-1)
    if TLnodelist(i,3)==TLnodelist(i+1,3) && TLnodelist(i,4)==TLnodelist(i+1,4)
        cnt = cnt + 1;
    else
        if cnt == 3
            bool_groupof3(i-2:i) = 1;
        end
        cnt = 1;
    end
end
            % the last else wasn't excuted
if cnt == 3
    bool_groupof3(i+1-2:i+1) = 1;
end
bool_groupof3 = logical(bool_groupof3);

tmp = TLnodelist(:,1); 
TLtri_group = tmp(bool_groupof3);
    
% 4. Find dihedral angles, first for each triangle within one group, then loop through the groups
    % for i = 1:#groups
        % for j = 1:3
            % get normal for the other two triangle in group, N1 and N2
            % find angle between N1 and N2, adjust to within 180
        % end
    % end
numGroup = length(TLtri_group)/3;
DihedralAngle = zeros(length(TLtri_group),1);
for i = 1:numGroup
    first = TLtri_group((i-1)*3 + 1);
    second = TLtri_group((i-1)*3 + 2);
    third = TLtri_group((i-1)*3 + 3);
    
        % acos doesn't work well for small angles
        % use atan2d(norm(cross(u,v)),dot(u,v))
    DihedralAngle(3*(i-1)+1) = atan2d(norm(cross(normal_full(second,:),normal_full(third,:))),dot(normal_full(second,:),normal_full(third,:)));
    DihedralAngle(3*(i-1)+2) = atan2d(norm(cross(normal_full(first,:),normal_full(third,:))),dot(normal_full(first,:),normal_full(third,:)));
    DihedralAngle(3*(i-1)+3) = atan2d(norm(cross(normal_full(first,:),normal_full(second,:))),dot(normal_full(first,:),normal_full(second,:)));
    
end
    

% Checkings
    % Within one group, check if angle between <n1,n2>, <n2,n3>, <n1,n3> sum to 360
sumGroupAngle = zeros(length(numGroup),1);
for i = 1:numGroup
    sumGroupAngle(i) = sum(DihedralAngle(((i-1)*3+1):((i-1)*3+3)));
end
    % check if n1, n2, n3 are coplanar
i = 923;

first = TLtri_group((i-1)*3 + 1);
second = TLtri_group((i-1)*3 + 2);
third = TLtri_group((i-1)*3 + 3);

n1 = normal_full(first,:);
n2 = normal_full(second,:);
n3 = normal_full(third,:);

tmp = cross(n1,n2);
dotproduct = dot(tmp,n3);


% write to txt file for D3D. 
fileID = fopen('060216_DAlist_Ni.txt','w');
fprintf(fileID,'%8d\n',length(DihedralAngle));
for i = 1 : length(DihedralAngle)
    %     NOTICE ID-1!
    fprintf(fileID,'%8d    %6.2f\n', (TLtri_group(i)-1) ,DihedralAngle(i));
end
fclose(fileID);



    
    