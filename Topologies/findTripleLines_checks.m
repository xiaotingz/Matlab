%% ##### Visualization #####
xdat = double(h5read(file, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'));
tri_not3group = triNodes(TLnodes(~mask_groupOf3,1),:);
FLs_not3group = FLs(TLnodes(~mask_groupOf3,1),:);
objGs = 1:4;
figure;
trisurf( tri_not3group(objGs,:), xdat(1,:), xdat(2,:), xdat(3,:) );
axis equal; axis off; 
rotate3d on

FL_group = FLs_not3group(objGs,:);
group_grains = unique(FL_group);
surf_grain = cell(length(group_grains),1);
for i = 1:length(group_grains)
    mask_grain = any(FLs_raw == group_grains(i), 2);
    surf_grain{i} = triNodes_raw(mask_grain, :);
end
trisurf(surf_grain{1}, xdat(1,:), xdat(2,:), xdat(3,:) );

% -----------------------------------------------------------------
% - What are the triple nodes/edges being shared by the 4 group triangles?
% - What exactly are the possible classes of 4-group triangles? In other
% words, all the cases of 2/4 grain groups, are how many are shared by 2
% grains and how many are shared by 4 grains? 
% -----------------------------------------------------------------


%% ##### Check #####
% - What's the triangle groups that didn't have 3 members? 
% - Are they all in groups of 2/4? Yes
% - How many 4-neighbor TLs? 
% - Except the manual-duplication induced false triple lines, how many 2-neighbor TLs?? 

% ##### non-3-memeber Triangle Groups #####
% """
%     - non3Groups = [ID_groupStart_TLnodes, groupSize]
% """
idx = 1;
cnt = 0;
non3Groups = [];
for i = 1:(length(mask_groupOf3) - 1)
    if mask_groupOf3(i+1) == 0
        cnt = cnt + 1;
    elseif mask_groupOf3(i) == 0 && mask_groupOf3(i+1) == 1
        non3Groups(idx, 1) = i - cnt + 1;
        non3Groups(idx, 2) = cnt;
        cnt = 0;
        idx = idx + 1;       
    end
end
if mask_groupOf3(i+1) == 0
    non3Groups(idx, 1) = i - cnt + 1;
    non3Groups(idx, 2) = cnt;
end

% ##### break large groups into elementary components #####
% --- Get the largeGroups ---
largeGroups = non3Groups((non3Groups(:,2) ~= 2 & non3Groups(:,2) ~= 4), :);

% --- Decompose the large groups ---
% """
%     - LGroups_decompose = [ID_groupStart_TLnodes, elementSize, groupSize]
% """
LGoups_decompose = [];
idx = 1;
cnt = 1;
for i = 1:length(largeGroups)
    group = TLnodes(largeGroups(i,1) : (largeGroups(i,1)+largeGroups(i,2)-1), :);
    for j = 1:length(group)-1
        if group(j,2) == group(j+1,2) && group(j,3) == group(j+1,3)
            cnt = cnt + 1;
        else
            LGoups_decompose(idx,1) = largeGroups(i,1) + j - cnt;
            LGoups_decompose(idx,2) = cnt;
            LGoups_decompose(idx,3) = largeGroups(i,2);
            cnt = 1;
            idx = idx + 1;
        end
    end
    LGoups_decompose(idx,1) = largeGroups(i,1) + (j + 1) - cnt;
    LGoups_decompose(idx,2) = cnt;
    LGoups_decompose(idx,3) = largeGroups(i,2);
    cnt = 1;
    idx = idx + 1;
end
if sum(LGoups_decompose(:,2)) ~= sum(largeGroups(:,2))
    warning('Eecompose of the large groups are wrong!');
end

% ##### quadLine fraction #####
numQL_tri = sum(non3Groups(:,2)==4) + sum(LGoups_decompose(:,2)==4);

% ##### quadLine length #####
% --- 1. If quadLines can have length>1: two QLs sharing one node ---
longQL = {};
cnt = 1;
for i = 1 : length(LGoups_decompose) - 1
    % --- if two 4-member-element sits in the same group & near each other ---
    % --- > 2 QLs will be recorded by segments of 2 ---
    if LGoups_decompose(i,2) == 4 && LGoups_decompose(i+1,2) == 4 && (LGoups_decompose(i,1) == (LGoups_decompose(i+1,1) - 4))
        % --- then check the nodes for connectivity --- 
        QLnodes = TLnodes(LGoups_decompose(i,1):LGoups_decompose(i,1)+7, 2:3);
        if length(unique(QLnodes)) < 4
            warning(['QLlength > 1 at:  ', num2str(LGoups_decompose(i,1))]);
            longQL{cnt, 1} = LGoups_decompose(i,1);
            idx_longQL_FL = TLnodes(LGoups_decompose(i,1):LGoups_decompose(i,1)+3, 1);
            longQL{cnt, 2} = unique(unique(FLs(idx_longQL_FL, :)));
            idx_longQL_FL = TLnodes(LGoups_decompose(i,1)+4:LGoups_decompose(i,1)+7, 1);
            longQL{cnt, 3} = unique(unique(FLs(idx_longQL_FL, :)));
            cnt = cnt + 1;
        end
    end
end
% --- Visualization ---
% idx = 1;
% tmp = TLnodes(longQL{idx,1}:longQL{idx,1}+7, :);
% idx_tmp = tmp(:,1);
% TLnodes_raw(idx_tmp, :);
% FLs(idx_tmp,:);
% objTris = triNodes(idx_tmp,:);

% trisurf( objTris, xdat(1,:), xdat(2,:), xdat(3,:) );
% rotate3d on

%%
% ##### number fraction Checks #####
% --- from the point of triangles ---
numNon3G_tri = length(mask_groupOf3) - sum(mask_groupOf3);
    % --- because I triplicated the corner triangles, introduced extra 2-member groups ---
numNon3G_tri_D3D = numNon3G_tri - (length(TLnodes) - length(TLnodes_raw))/2;
numQL_tri = (sum(non3Groups(:,2)==4) + sum(LGoups_decompose(:,2)==4))*4;
numTwoL_tri = (sum(non3Groups(:,2)==2) + sum(LGoups_decompose(:,2)==2))*2 - ...
    (length(TLnodes) - length(TLnodes_raw))/2;
if numTwoL_tri + numQL_tri ~= numNon3G_tri_D3D
    warning('Hey the non-3 group triangles numbers are wrong!');
end

% --- from the point of segments ---
numTLs = sum(mask_groupOf3)/3
sumQLs = sum(non3Groups(:,2)==4) + sum(LGoups_decompose(:,2)==4)
    % --- numTriplicatedTris = (length(TLnodes) - length(TLnodes_raw))/2 ---
    % --- theInducedGroup = numTriplicatedTris/2 ---
sumTwoLs =  sum(non3Groups(:,2)==2) + sum(LGoups_decompose(:,2)==2) - (length(TLnodes) - length(TLnodes_raw))/4

