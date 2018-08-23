function triple_line = findTripleLines(file)
% ############################################################################
% - tripleLines = [n, 3], n is the number of triple lines in the volume, 3
%       stands for the threeGrainIDs encapsulating the triple line.
% ############################################################################
% ----------------------- load debug data -----------------------
% clear
% file = '/Users/xiaotingzhong/Desktop/Datas/synthetic/180502_CubicSingleEquiaxedOut.dream3d';
% file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_mesh_mergeTwin.dream3d';
% ---------------------------------------------------------------

% ##### Read in data ##### 
tri_node_raw = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
face_lable_raw = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
node_type = double(h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType'))';

% ##### Find the triangles on tripleline #####
%     """
%     - treat QuadNodes as TLnodes
%     - the triple line triangles are those which have at least 2 TLnodes or QuadNodes 
%            TLnodes=3 happens for the triangles sitting on a corner 
%     """
node_type(node_type==4) = 3;
mask_TLtri = (sum(node_type(tri_node_raw)==3, 2) >= 2);
FLs = face_lable_raw(mask_TLtri, :);
triNodes = tri_node_raw(mask_TLtri, :);

    
% ##### Sort the triangles into groups of 3 #####
%     """
%     - creat a TLnode list to record the triple line nodes of each triangle,
%           group triangles according to the common triple line node 
%     - special case are the triangles having 3 quad nodes, then it's difficult to tell which are the two nodes in use
%           in such case, duplicate the triangle
%     """
mask_TLnodes = (node_type(triNodes)==3);
TLnodes_raw = triNodes.*mask_TLnodes;
TLnodes_raw = sort(TLnodes_raw,2);
% TLnodes(:,1) = (1:length(TLnodes))';
mask_3quadTris = ~any(TLnodes_raw == 0, 2);
TLnodes = zeros(length(TLnodes_raw) + sum(mask_3quadTris)*2, 3);
idx = 1;
for i = 1:length(TLnodes_raw)
    if TLnodes_raw(i, 1) == 0
        TLnodes(idx, 1) = i;
        TLnodes(idx, 2:3) = TLnodes_raw(i, 2:3);
        idx = idx + 1;
    else
        TLnodes(idx:idx+2, 1) = i;
        TLnodes(idx, 2:3) = TLnodes_raw(i, 1:2);
        TLnodes(idx+1, 2:3) = TLnodes_raw(i, 2:3);
        TLnodes(idx+2, 2) = TLnodes_raw(i, 1);
        TLnodes(idx+2, 3) = TLnodes_raw(i, 3);
        idx = idx + 3;
    end
end
% ---- TLnodes = [IDinTLnodes_raw, TLnode1, TLnode2] ----
TLnodes = sortrows(TLnodes,[2,3]);

% --- if a group contain more or less than 3 triangles, ignore ---
mask_groupOf3 = zeros(length(TLnodes),1);
cnt = 1;
for i = 1:(length(TLnodes)-1)
    if TLnodes(i,2)==TLnodes(i+1,2) && TLnodes(i,3)==TLnodes(i+1,3)
        cnt = cnt + 1;
    else
        if cnt == 3
            mask_groupOf3(i-2:i) = 1;
        end
        cnt = 1;
    end
end
% --- the last else wasn't excuted ---
if cnt == 3
    mask_groupOf3(i+1-2:i+1) = 1;
end
mask_groupOf3 = logical(mask_groupOf3);

% ##### get triple lines from the triangle groups #####
TLtri_FLs= FLs(TLnodes(mask_groupOf3, 1), :);
if rem(length(TLtri_FLs), 3) ~= 0 
    warning('ERROR, length(TLtri_FLs)/3 is not integer')
end
triple_line = zeros(length(TLtri_FLs)/3, 3);
for i = 1:length(TLtri_FLs)/3
    FLgroup = TLtri_FLs((i-1)*3+1 : (i-1)*3+3, :);
    if length(unique(FLgroup)) == 3
        triple_line(i,:) = unique(FLgroup);
    else
        warning(['length(unique(FLgroup)) ~= 3 at FLgroup=', num2str(i)]);
    end
end
triple_line = sort(triple_line, 2);
triple_line = unique(triple_line, 'rows');

% --- clean the data --- 
mask_surfTL = any(triple_line==0, 2);
triple_line = triple_line(~mask_surfTL, :);

end






