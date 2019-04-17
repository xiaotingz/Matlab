% function dist_ctwin = calcDistFromCTwin(file, fl_obj, O)
% ##########################################################################
% * Input
%     - fls = [n , 2]
%           Labels of the faces, for which misorientation angle needs to be calcualted.
%     - O = [3, 3, 24]
%           The symmetry operators.
% * Output
%     - dist_twin = [n, 1]
%           A weighted distance between given grain boundary (five
%           parameters) and coherent twin (60@[111], (111)) 
% * Note
%     - For the ith grain face, check all triangles on it and for every
%     tri_(i,j), calculate its distance to the coherent twin by Morawiec's
%     metric and weight these distances by triangle area. 
%     - Morawiec's formulas in FIND_NEIGH.f90, using n_crys = n.T*g
% ##########################################################################
% ----------------------- load debug data -----------------------
file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin2_smooth.dream3d');
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/181107.mat', 'tracked_uniqueface_an4')
fl_obj = tracked_uniqueface_an4;
clear tracked_uniqueface_an4
load('/Users/xiaotingzhong/Documents/Matlab/D3Ddistris/crysSym.mat');
% ---------------------------------------------------------------

% ######################### Prepare Data #########################
% ----- load data -----
ea =  h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEAs').';
fl =  h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels').';
n = h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals').';
area = h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas').';
ea(1, :) = [];
ea = rad2deg(ea);
mask = all(fl > 0, 2);
fl = fl(mask, :);
n = n(mask, :);
area = area(mask);

% ----- unify face label -----
fl_obj = sort(fl_obj, 2);
fl = sort(fl, 2);

% ----- prepare data for twin -----
dg_twin = AAToG(60, [1, 1, 1]);
n_twin = [1, 1, 1]';
n_twin = n_twin / norm(n_twin);
b_twin = calcMorawiecBMat(dg_twin, n_twin);

% ----- prepare symmetry operator for morawiec matrix-----
bsym = zeros(4, 4, size(O, 3));
for i = 1:size(O, 3)
    bsym(1:3, 1:3, i) = O(:,:,i);
    bsym(4, 4, i) = 1;
end

% %%
% --------------------------------------------------------------------------------------
% ----- output to python -----
% idx_tri = 1;
% g_1 = EAtoG(ea(obj_fl(i, 1), :));
% g_2 = EAtoG(ea(obj_fl(i, 2), :));
% normal = n(idx_tri, :);
% 
% n_1 = (normal * g_1)';
% n_2 = (normal * g_2)';
% dg = g_1'*g_2;
% b_1 = calcMorawiecBMat(dg, n_1);
% b_2 = calcMorawiecBMat(dg', -n_2);
% spv4dist(b_1, b_2, bsym)
% 
% fileID = fopen('tmp.txt', 'w');
% fprintf(fileID, '%s%6.3f, %6.3f, %6.3f%s\n', 'g1 = np.array([[', g_1(1,:), '],');
% fprintf(fileID, '\t\t%s%6.3f, %6.3f, %6.3f%s\n', '[', g_1(2,:), '],');
% fprintf(fileID, '\t\t%s%6.3f, %6.3f, %6.3f%s\n', '[', g_1(3,:), ']])');
% fprintf(fileID, '%s%6.3f, %6.3f, %6.3f%s\n','g2 = np.array([[',  g_2(1,:), '],');
% fprintf(fileID, '\t\t%s%6.3f, %6.3f, %6.3f%s\n', '[', g_2(2,:), '],');
% fprintf(fileID, '\t\t%s%6.3f, %6.3f, %6.3f%s\n', '[', g_2(3,:), ']])');
% fprintf(fileID, '%s%6.3f, %6.3f, %6.3f%s\n', 'n = np.array([', normal, '])');
% fclose(fileID)
% --------------------------------------------------------------------------------------
% %%

%%
% #################### Calculate dist_ctwin For Each Grain Face Of Interest ####################
% """
% cn = n'*g
%     Usually cn = g*n, however, Morawiec is using a different formula for b also. 
%     Probably fine as long as consistent. 
%     See FIND_NEIGH.f90, SUBROUTINE find_nearest. (imain and jmain are for TJs, ib and jb are for three GBs at a TJ)
% """
dist_ctwin = zeros(length(fl_obj), 1);
for i = 1:length(fl_obj)
    % ----- Calculate misorientation across face -----
    g_1 = EAtoG(ea(fl_obj(i, 1), :));
    g_2 = EAtoG(ea(fl_obj(i, 2), :));
    dg = g_1 * g_2';

    % ----- Get triangles on the grain face -----
    mask = ismember(fl, fl_obj(i, :), 'rows');
    n_face = n(mask, :);
    area_face = area(mask, :);

    % ----- For each triangle, calculate distance to ctwin -----
    dist_face = zeros(size(area_face));
    for j = 1:sum(mask)
        % """
        % calcMorawiecBMat(dg, n_1) = calcMorawiecBMat(dg', n_2)
        % """
        n_1 = (n_face(j, :) * g_1)';
        b = calcMorawiecBMat(dg, n_1);
        
        tmp_dist_1 = spv4dist(b, b_twin, bsym);
        tmp_dist_2 = spv4dist(b_twin, b, bsym);
        
        dist_face(j) = min(tmp_dist_1, tmp_dist_2);        
    end
    
    % ----- Calculate weighted distance of the entire grain face -----
    dist_ctwin(i) = sum(dist_face .* area_face) / sum(area_face);
    
end



    







