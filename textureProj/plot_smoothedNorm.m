% 1. Read in data
tri_full = textread('/Users/xiaotingzhong/Desktop/Datas/Ferrite/Jan.31 Ferrite/Jan31_F_triangle.txt');
node_full = textread('/Users/xiaotingzhong/Desktop/Datas/Ferrite/Jan.31 Ferrite/Jan31_F_node.txt');
normal_full = roundn(h5read('/Users/xiaotingzhong/Desktop/Datas/Ferrite/Jan.31 Ferrite/Jan31_Fc_stats.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals'),-5);

    % adjust count starting from 0 to 1
node_full(1,:) = [];
node_full(:,1) = node_full(:,1) + 1;
tri_full(1,:) = [];
tri_full(:,1:4) = tri_full(:,1:4) + 1;
tri_full(:,5:7) = [];
normal_full = normal_full.';

% 2. Find the triangles on tripleline
    % first clean surface triangles, deal normal list at the same time
tmp1 = tri_full(:,1); tmp2 = tri_full(:,2); tmp3 = tri_full(:,3); 
tmp4 = tri_full(:,4); tmp5 = tri_full(:,5); tmp6 = tri_full(:,6);
bool1 = (tri_full(:,5) > 0 & tri_full(:,6) > 0);
GBtri(:,1) = tmp1(bool1); GBtri(:,2) = tmp2(bool1); GBtri(:,3) = tmp3(bool1);
GBtri(:,4) = tmp4(bool1); GBtri(:,5) = tmp5(bool1); GBtri(:,6) = tmp6(bool1);

tmp7 = normal_full(:,1); tmp8 = normal_full(:,2); tmp9 = normal_full(:,3); 
normals(:,1) = tmp7(bool1); normals(:,2) = tmp8(bool1); normals(:,3) = tmp9(bool1);
clear tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9 bool1

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
tmp7 = normals(:,1); tmp8 = normals(:,2); tmp9 = normals(:,3); 
clear normals
normals(:,1) = tmp7(bool2); normals(:,2) = tmp8(bool2); normals(:,3) = tmp9(bool2); 
clear tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9 bool2 
    
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

load('smoothed_Fnorm_050316.mat');
normals = sortrows(normals);

TLnodelist = sort(TLnodelist,2);
TLnodelist = [TLtri(:,1),zeros(length(TLtri),1),TLnodelist,zeros(length(TLtri),1),normals];
TLnodelist = sortrows(TLnodelist,3);

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

tmp1 = TLnodelist(:,7); tmp2 = TLnodelist(:,8); tmp3 = TLnodelist(:,9);
groupedNorm(:,1) = tmp1(bool_groupof3); 
groupedNorm(:,2) = tmp2(bool_groupof3); 
groupedNorm(:,3) = tmp3(bool_groupof3);

numGroup = length(groupedNorm)/3;
DA = zeros(length(groupedNorm),1);
for i = 1:numGroup
    first = (i-1)*3 + 1;
    second = (i-1)*3 + 2;
    third = (i-1)*3 + 3;
    
        % acos doesn't work well for small angles
        % use atan2(norm(cross(u,v)),dot(u,v))
    DA(first) = atan2d(norm(cross(groupedNorm(second,:),groupedNorm(third,:))),dot(groupedNorm(second,:),groupedNorm(third,:)));
    DA(second) = atan2d(norm(cross(groupedNorm(first,:),groupedNorm(third,:))),dot(groupedNorm(first,:),groupedNorm(third,:)));
    DA(third) = atan2d(norm(cross(groupedNorm(first,:),groupedNorm(second,:))),dot(groupedNorm(first,:),groupedNorm(second,:)));
  
end
    







