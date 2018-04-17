% 1. Read in data
% traingle_full(triangle_ID, coor1_ID, coor2_ID, coor2_ID, --, --, --, label_1, label_2)
tri_full = textread('/Users/xiaotingzhong/Desktop/Datas/Jan.31 Austenite/Jan31_A_croptriangle.txt');
node_full = textread('/Users/xiaotingzhong/Desktop/Datas/Jan.31 Austenite/Jan31_A_cropnode.txt');

    % adjust count starting from 0 to 1
node_full(1,:) = [];
node_full(:,1) = node_full(:,1) + 1;
tri_full(1,:) = [];
tri_full(:,1:4) = tri_full(:,1:4) + 1;
tri_full(:,5:7) = [];

% 2. Find the triangles on tripleline
    % first clean surface triangles
bool1 = (tri_full(:,5) > 0 & tri_full(:,6) > 0);
Innertri = tri_full(bool1, :);

    % then check node type, if 2 of the three nodes are tripleline node (type3), keep
bool2 = zeros(length(Innertri),1);
for i = 1:length(Innertri)
    cnt = 0;
    for j = 2:4
        nodeType = node_full(Innertri(i,j),2);
        if nodeType == 3
            cnt = cnt + 1;
        end
    end
    
    if cnt == 2
        bool2(i) = 1;
    end
end
bool2 = logical(bool2);
TLtri = Innertri(bool2, :);
    
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
        % TLnodelist = [triangle_ID, 0, triangleTripleLineNode1, triangleTripleLineNode2]        
    % first sort the smaller label to be at front thus unify the labels for
    % triangles on the same grain face. 
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
% In the data TLtri_group, IDs of the three triangles sharing one triple
% line is organized together. From the triangle_ID should be able to get
% the face_label then face_ID 
TLtri_group = tmp(bool_groupof3);