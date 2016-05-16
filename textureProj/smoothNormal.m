% 1. Read in data
% % Austenite
% tri_full = textread('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_A_triangle.txt');
% node_full = textread('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_A_node.txt');
% normal_full = roundn(h5read('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_A_stats.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals'),-5);
% % Ferrite
tri_full = textread('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_F_triangle.txt');
node_full = textread('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_F_node.txt');
normal_full = roundn(h5read('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_Fc_stats.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals'),-5);
% % Cropped
% tri_full = textread('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_A_croptriangle.txt');
% node_full = textread('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_A_cropnode.txt');
% normal_full = roundn(h5read('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_A_cropstats.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals'),-5);

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
GoodTri(:,1) = tmp1(bool1); GoodTri(:,2) = tmp2(bool1); GoodTri(:,3) = tmp3(bool1);
GoodTri(:,4) = tmp4(bool1); GoodTri(:,5) = tmp5(bool1); GoodTri(:,6) = tmp6(bool1);
clear tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 bool1

    % then check node type, if 2 of the three nodes are tripleline node (type3), keep
bool2 = zeros(length(GoodTri),1);
for i = 1:length(GoodTri)
    cnt = 0;
    for j = 2:4
        nodeType = node_full(GoodTri(i,j),2);
        if nodeType == 3
            cnt = cnt + 1;
        end
    end
    
    if cnt == 2
        bool2(i) = 1;
    end
end
bool2 = logical(bool2);
tmp1 = GoodTri(:,1); tmp2 = GoodTri(:,2); tmp3 = GoodTri(:,3); 
tmp4 = GoodTri(:,4); tmp5 = GoodTri(:,5); tmp6 = GoodTri(:,6);
TLtri(:,1) = tmp1(bool2); TLtri(:,2) = tmp2(bool2); TLtri(:,3) = tmp3(bool2);
TLtri(:,4) = tmp4(bool2); TLtri(:,5) = tmp5(bool2); TLtri(:,6) = tmp6(bool2);
clear tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 bool2 


% 3. Smooth triangle normal with its nearest neighbor
    % sort the triangles according to faceLabel list
GoodTri(:,5:6) = sort(GoodTri(:,5:6),2);
GoodTri = sortrows(GoodTri,[5 6]); 
GoodTri = [GoodTri(:,1),zeros(length(GoodTri),1),GoodTri(:,2:6)];

FaceTri_TL = sort(TLtri(:,5:6),2);
FaceTri_TL = [TLtri(:,1),ones(length(FaceTri_TL),1),FaceTri_TL];
FaceTri_TL = sortrows(FaceTri_TL,[3 4]);

    % identify triple line triangles among the face triangles
cnt = 1;
for i = 1:length(GoodTri)
    if GoodTri(i) == FaceTri_TL(cnt)
        GoodTri(i,2) = 1;
        cnt = cnt + 1;
    end
end

    % count #triangles in each face
            % count #triangles for only faces with triple line
        % for i = 1:length(GoodTri)
        %     if GoodTri(i,6) == Faces(cnt_f,1) && GoodTri(i,7) == Faces(cnt_f,2)
        %         cnt_t = cnt_t + 1;
        %     elseif GoodTri(i+1,6) == Faces(cnt_f+1,1) && GoodTri(i+1,7) == Faces(cnt_f+1,2)
        %         Faces(cnt_f,3) = cnt_t;
        %         cnt_f = cnt_f + 1;
        %         cnt_t = 1 + adjust;
        %         adjust = 0;
        %     else 
        %         adjust = -1;
        %         continue
        %     end
        % end
        % Faces(cnt_f,3) = cnt_t;
Faces = unique(GoodTri(:,6:7),'rows');
Faces = [Faces,zeros(length(Faces),1)];
cnt_f = 1;
cnt_t = 0;
for i = 1:length(GoodTri)
    if GoodTri(i,6) == Faces(cnt_f,1) && GoodTri(i,7) == Faces(cnt_f,2)
        cnt_t = cnt_t + 1;
    else
        Faces(cnt_f,3) = cnt_t;
        cnt_f = cnt_f + 1;
        cnt_t = 1;
    end
end
Faces(cnt_f,3) = cnt_t;
clear cnt_f cnt_t
    
    % for i = 1:#Faces 
        % for j = 1:#(triangle in the ith Faces)
            % check if is a triple line triangle
                % for k = 1:#(triangle in the ith Faces)
                    % check common vertex, check if nearest neighbors (NN_j)
                    % if NN_jcheck original label, adjust normal direction and record
                % end
                % average normal of the (k+1) triangles in the patch
            % end
        % end
    % end   
cnt_patch = 1;
cnt_tl = 1;
begin = 1;
final = 0;
normals = zeros(length(TLtri),4);
for i = 1:length(Faces)
    begin = 1 + sum(Faces(1:(i-1),3));
    final = sum(Faces(1:i,3));
    for j = begin:final
        if GoodTri(j,2) == 1
            vertex_j = GoodTri(j,3:5);
            for k = begin:final
                if ~isempty(intersect(vertex_j,GoodTri(k,3:5)))
                    patch_norm_j(cnt_patch,:) = normal_full(GoodTri(k,1),:);
                    if tri_full(GoodTri(k,1),5:6) ~= tri_full(GoodTri(j,1),5:6)
                        patch_norm_j(cnt_patch,:) = -patch_norm_j(cnt_patch);
                    end
                    cnt_patch = cnt_patch + 1;
                end
            end
            patch_aveNorm = sum(patch_norm_j);
            patch_aveNorm = patch_aveNorm/norm(patch_aveNorm);
            normals(cnt_tl,1) = GoodTri(j,1);
            normals(cnt_tl,2:4) = patch_aveNorm;
            cnt_tl = cnt_tl + 1;
            cnt_patch = 1;
        end
    end
end
            
            
 
% 4. Sort the triangles into groups of 3
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
    % if the group contain more than 3 triangles, ignore
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


% 5. Find Dihedral Angles
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
        % use atan2(norm(cross(u,v)),dot(u,v))
    DihedralAngle(3*(i-1)+1) = atan2d(norm(cross(normal_full(second,:),normal_full(third,:))),dot(normal_full(second,:),normal_full(third,:)));
    DihedralAngle(3*(i-1)+2) = atan2d(norm(cross(normal_full(first,:),normal_full(third,:))),dot(normal_full(first,:),normal_full(third,:)));
    DihedralAngle(3*(i-1)+3) = atan2d(norm(cross(normal_full(first,:),normal_full(second,:))),dot(normal_full(first,:),normal_full(second,:)));
    
end

    
    
    
    