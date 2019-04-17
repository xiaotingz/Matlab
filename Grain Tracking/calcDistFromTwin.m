% function dist_twin = calcDistFromTwin(file, fl_obj, O)
% ##########################################################################
% * Input
%     - fls = [n , 2]
%           labels of the faces, for which misorientation angle needs to be calcualted
%     - O = [3, 3, 24]
%           symmetry operators
% * Output
%     - dist_twin = [n, 1]
%           distance between given misorientations and sigma3.
% * Note
%     - This file calculates the misorientation angle across a given boundary.
%     - In theory, [A, B] and [B, A] should have the same misorientation
%     angle. In practice, there are small numeric problems but the
%     deviation is really small.
% ##########################################################################
% ----------------------- load debug data -----------------------
file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin2_smooth.dream3d');
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/181107.mat', 'tracked_uniqueface_an4')
fl_obj = tracked_uniqueface_an4;
clear tracked_uniqueface_an4
load('/Users/xiaotingzhong/Documents/Matlab/D3Ddistris/crysSym.mat');
% ---------------------------------------------------------------

% ############################### Prepare Data ###############################
ea =  h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEAs').';
ea(1, :) = [];
ea = rad2deg(ea);

% ----- dg_twin in FZ -----
dg_twin = AAToG(60, [1, 1, 1]);

%%
% ############################### Calc MisA(dg, dg_twin) ###############################
% """
% - Distance between misorientations is calculated as 1) put the two misorientations in FZ.
% 2) calc misorientation angle between two misorientations. 
% - misA formula: Kryz' thesis, P21.
% - Similar calculation: dgInFZ.m
% """
% ----- initialize the result -----
dist_twin = ones(size(fl_obj, 1), 1)*180;

% for i = 1:size(fl_obj, 1)
% for i = randi(size(fl_obj, 1), 1)
%     disp(['face: ', int2str(i)])
%     disp(mis_ang_obj(i));
%     
%     g_1 = EAtoG(ea(fl_obj(i, 1), :));
%     g_2 = EAtoG(ea(fl_obj(i, 2), :));
%     min_ang = 180;
%     for j = 1:24
%         for k = 1:24
%             gg_1 = O(:,:,j)*g_1;
%             gg_2 = O(:,:,k)*g_2;
%             
%             dg = gg_1*gg_2';
% %             if 1 - abs(0.5*(trace(dg)-1)) > 0.001
%                 ang = acosd(0.5*(trace(dg)-1));
%                 if ang < min_ang
%                     min_ang = ang;
%                 end
% %             else
% %                 disp(['exit1, at j,k = ', num2str([j,k])])
% %                 min_ang = 0.0;
% %             end
%                 
%             
% %             if 1 - abs(0.5*(trace(dg)-1)) > 0.001
%                 dg = gg_2*gg_1';
%                 ang = acosd(0.5*(trace(dg)-1));
%                 if ang < min_ang
%                     min_ang = ang;
%                 end
% %             else
% %                 disp(['exit2, at j,k = ', num2str([j,k])])
% %                 min_ang = 0.0;
% %             end
%         end
%     end
%     
%     disp(min_ang)
%     

for i = tmp(533)
    g_1 = EAtoG(ea(fl_obj(i, 1), :));
    g_2 = EAtoG(ea(fl_obj(i, 2), :));

    % ----- First put the misorientation in FZ -----
    [~, ~, dg_fz] = dgInFZ(g_1, g_2, O);
    
    % ----- Then treat misorientations as orientations to calculate misorientations angle (of misorientation) -----
    dg_1 = dg_fz;
    dg_2 = dg_twin;
    for j = 1:24
        for k = 1:24
            gdg_1 = O(:,:,j)*dg_1;
            gdg_2 = O(:,:,k)*dg_2;
            
            ddg = gdg_1*gdg_2';
            ang = acosd(0.5*(trace(ddg)-1));
            if ang < dist_twin(i)
                dist_twin(i) = ang;
            end
                
            ddg = gdg_2*gdg_1';
            ang = acosd(0.5*(trace(ddg)-1));
            if ang < dist_twin(i)
                dist_twin(i) = ang;
            end
        end
    end
    disp([mis_ang_obj(i), dist_twin(i)])
end











%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% """
% Check the calculation with misorientation angle of the boundary itself. 
% """
% ##### Load data #####
num_neigh =  h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors').';
neigh_list =  h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList');
mis_ang = h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/MisorientationList');
num_neigh(1) = [];


% ##### Prepare index #####
fl_full = zeros(size(neigh_list));
for i = 1:size(num_neigh)
    idx_s = sum(num_neigh(1 : i-1)) + 1;
    idx_e = sum(num_neigh(1 : i));
    fl_full(idx_s : idx_e) = i;
end
fl_full = [fl_full, neigh_list];

% ##### Get misorientation distance, as misorientation angle, to twin #####
mask = ismember(fl_full, fl_obj, 'row');
mis_ang_tmp = mis_ang(mask);
idx_helper = (1:size(fl_obj, 1))';
fl_idx_tmp = [fl_obj, idx_helper];
fl_idx_tmp = sortrows(fl_idx_tmp, [1, 2]);
if sum(all(fl_idx_tmp(:, 1:2) == fl_full(mask, :), 2)) ~= size(fl_obj, 1)
    warning('calcDistFromTwin.m, index not right!')
end
idx_convert = sortrows([idx_helper, fl_idx_tmp(:, 3)], 2);
mis_ang_obj = mis_ang_tmp(idx_convert(:,1));


%%
% dist_twin = dist_twin_calc;
mask_misA60 = abs(mis_ang_obj - 60) < 5;
mask_closetwin = abs(dist_twin) < 5;

disp([sum(mask_misA60), sum(mask_closetwin), sum(mask_misA60==mask_closetwin & mask_closetwin==1)]);

mask_check = (mask_closetwin & mask_misA60==0);
tmp = [idx_helper(mask_check), mis_ang_obj(mask_check), dist_twin(mask_check)];








