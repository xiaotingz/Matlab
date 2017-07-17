% load data
clear
file = ('/Users/xiaotingzhong/Desktop/091616_STO_1470/101116_V4_misA3/sub1_misA3_recons_curv.dream3d');

num_of_neigh = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
neighborList = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NeighborList'));
DisAngle = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/MisorientationList'));
facelabel = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels')).';
    % notice use the ABS(triCurv)
triCurv_raw = abs(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));
triangle_area_raw = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-8);
surf_grain = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/SurfaceFields'));

% get the disA of the unique faces 
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

face_disA_partial = [A, neighborList, DisAngle];
for i = 1:length(face_disA_partial)
    if face_disA_partial(i,1) > face_disA_partial(i,2)
        tmp = face_disA_partial(i,1);
        face_disA_partial(i,1) = face_disA_partial(i,2);
        face_disA_partial(i,2) = tmp;
    end
end
face_disA_partial = sortrows(face_disA_partial, [1,2]);

face_disA = zeros(length(neighborList)/2,3);
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

% sort the triangles in the face way
bool1 = facelabel(:,1) > 0 & facelabel(:,2) > 0 & ~isnan(triCurv_raw(:)) & triCurv_raw(:) < 100 & triCurv_raw(:) > -100;
tmp1 = facelabel(:,1);
tmp2 = facelabel(:,2);
triInfo(:,1) = tmp1(bool1);
triInfo(:,2) = tmp2(bool1);
triInfo(:,3) = triangle_area_raw(bool1);
triInfo(:,4) = triCurv_raw(bool1);
clear tmp1 tmp2

for i = 1:length(triInfo)
    if triInfo(i,1) > triInfo(i,2)
        tmp = triInfo(i,1);
        triInfo(i,1) = triInfo(i,2);
        triInfo(i,2) = tmp;
    end
end
triInfo = sortrows(triInfo,[1,2]);

face_curv = zeros(length(neighborList)/2,4);
ca_product = 0;
area = 0;
cnt = 1;
for i = 1:length(triInfo)-1
    if triInfo(i,1) == triInfo(i+1,1) && triInfo(i,2) == triInfo(i+1,2)
        ca_product = ca_product + triInfo(i,3)*triInfo(i,4);
        area = area + triInfo(i,3);
    else
        ca_product = ca_product + triInfo(i,3)*triInfo(i,4);
        area = area + triInfo(i,3);
        
        face_curv(cnt,1) = triInfo(i,1);
        face_curv(cnt,2) = triInfo(i,2);
        face_curv(cnt,3) = ca_product;
        face_curv(cnt,4) = area;
        ca_product = 0;
        area = 0;
        cnt = cnt + 1;
    end
end
ca_product = ca_product + triInfo(i+1,3)*triInfo(i+1,4);
area = area + triInfo(i+1,3);
face_curv(cnt,1) = triInfo(i+1,1);
face_curv(cnt,2) = triInfo(i+1,2);
face_curv(cnt,3) = ca_product;
face_curv(cnt,4) = area;
        
%% average curv according to area, including surface grains
start = 0;
nstep = 35;
stepsize = 2;
data_grid = zeros(nstep,5);
total_ca_product = 0;
total_area = 0;
cnt_faces = 0;
for j = 1 : nstep
    for i = 1 : length(face_disA)
        % data_final = [#triangles on the face, total curvature on the face, disorientation]
        if face_disA(i,3) > start && face_disA(i,3) <= start + stepsize
            total_ca_product = total_ca_product + face_curv(i,3);
            total_area = total_area + face_curv(i,4);
            cnt_faces = cnt_faces + 1;
        end
    end
    data_grid(j,1) = start + stepsize/2;
    data_grid(j,2) = total_ca_product;
    data_grid(j,3) = total_area;
    data_grid(j,4) = total_ca_product/total_area;
    data_grid(j,5) = cnt_faces;
    
    start = start + stepsize;
    total_ca_product = 0;
    total_area = 0;
    cnt_faces = 0;
end

%% average curv according to area, without surface grains
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

tmp3 = face_disA(:,3);
tmp4 = face_curv(:,3);
tmp5 = face_curv(:,4);
data_inbulk(:,1) = tmp3(bool_inbulk);
data_inbulk(:,2) = tmp4(bool_inbulk);
data_inbulk(:,3) = tmp5(bool_inbulk);
clear tmp3 tmp4 tmp5

start = 0;
nstep = 35;
stepsize = 2;
data_grid = zeros(nstep,5);
total_ca_product = 0;
total_area = 0;
cnt_faces = 0;
for j = 1 : nstep
    for i = 1 : length(data_inbulk)
        % data_final = [#triangles on the face, total curvature on the face, disorientation]
        if data_inbulk(i,1) > start && data_inbulk(i,1) <= start + stepsize
            total_ca_product = total_ca_product + data_inbulk(i,2);
            total_area = total_area + data_inbulk(i,3);
            cnt_faces = cnt_faces + 1;
        end
    end
    data_grid(j,1) = start + stepsize/2;
    data_grid(j,2) = total_ca_product;
    data_grid(j,3) = total_area;
    data_grid(j,4) = total_ca_product/total_area;
    data_grid(j,5) = cnt_faces;
    
    start = start + stepsize;
    total_ca_product = 0;
    total_area = 0;
    cnt_faces = 0;
end


%% average curv according to #faces, including surface grains
% start = 0;
% nstep = 35;
% stepsize = 2;
% data_grid = zeros(nstep,4);
% total_ca_product = 0;
% cnt_faces = 0;
% for j = 1 : nstep
%     for i = 1 : length(face_disA)
%         % data_final = [#triangles on the face, total curvature on the face, disorientation]
%         if face_disA(i,3) > start && face_disA(i,3) <= start + stepsize
%             total_ca_product = total_ca_product + face_curv(i,3)/face_curv(i,4);
%             cnt_faces = cnt_faces + 1;
%         end
%     end
%     data_grid(j,1) = start + stepsize/2;
%     data_grid(j,2) = total_ca_product;
%     data_grid(j,3) = cnt_faces;
%     data_grid(j,4) = total_ca_product/cnt_faces;
%     start = start + stepsize;
%     total_ca_product = 0;
%     cnt_faces = 0;
% end
% 
% 
% %% average curv according to #faces, without surface grains
% bool_inbulk = ones(length(face_disA),1);
% for i = 1:length(face_disA)
% %     the indexes in surface_grains are grainID+1
%     gA = face_disA(i,1) + 1;
%     gB = face_disA(i,2) + 1;
%     if surf_grain(gA) == 1 && surf_grain(gB) == 1
%         bool_inbulk(i) = 0;
%     end
% end
% bool_inbulk = logical(bool_inbulk);
% 
% tmp3 = face_disA(:,3);
% tmp4 = face_curv(:,3);
% tmp5 = face_curv(:,4);
% data_inbulk(:,1) = tmp3(bool_inbulk);
% data_inbulk(:,2) = tmp4(bool_inbulk);
% data_inbulk(:,3) = tmp5(bool_inbulk);
% clear tmp3 tmp4 tmp5
% 
% start = 0;
% nstep = 35;
% stepsize = 2;
% data_grid = zeros(nstep,4);
% total_ca_product = 0;
% cnt_faces = 0;
% for j = 1 : nstep
%     for i = 1 : length(data_inbulk)
%         % data_final = [#triangles on the face, total curvature on the face, disorientation]
%         if data_inbulk(i,1) > start && data_inbulk(i,1) <= start + stepsize
%             total_ca_product = total_ca_product + data_inbulk(i,2)/data_inbulk(i,3);
%             cnt_faces = cnt_faces + 1;
%         end
%     end
%     data_grid(j,1) = start + stepsize/2;
%     data_grid(j,2) = total_ca_product;
%     data_grid(j,3) = cnt_faces;
%     data_grid(j,4) = total_ca_product/cnt_faces;
%     start = start + stepsize;
%     total_ca_product = 0;
%     cnt_faces = 0;
% end        
%         
%         
%         
% 
%         
%         
%         
%         
%         