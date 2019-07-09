function dists = calcDistFromMisorientation(file, fl_obj, dg_obj)
% ##########################################################################
% * Input
%     - fls = [n , 2]
%           labels of the faces, for which misorientation angle needs to be calcualted
%     - dg_obj = [3, 3, n]
%           orientation matrix of the objective misorientations.
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
% load('/Volumes/XIAOTING/Ni/working/190425_Hsmooth_geo_topo_an5crop2.mat', ...
%     'tracked_uniqueface_an4', 'tracked_uniqueface_an5');
% file = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth.dream3d';
% file = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% fl_obj = tracked_uniqueface_an4_full;
% % fl_obj = tracked_uniqueface_an5;
% % clear tracked_uniqueface_an4
% 
% % ----- dg_twin in FZ -----
% dg_obj = zeros(3, 3, 6);
% dg_obj(:,:,1) = AAToG(60, [1, 1, 1]);     % sigma3
% dg_obj(:,:,2) = AAToG(36.86, [1, 0, 0]);  % sigma5
% dg_obj(:,:,3) = AAToG(38.21, [1, 1, 1]);  % sigma7
% dg_obj(:,:,4) = AAToG(38.94, [1, 1, 0]);  % sigma9
% dg_obj(:,:,5) = AAToG(31.59, [1, 1, 0]);  % sigma27a
% dg_obj(:,:,6) = AAToG(35.43, [2, 1, 0]);  % sigma27b
% ---------------------------------------------------------------

% ----- initialize the result variable -----
num_obj_miso = size(dg_obj, 3);
dists = ones(size(fl_obj, 1), num_obj_miso)*180;

% ############################### Prepare Data ###############################
ea =  h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles').';
ea(1, :) = [];
ea = rad2deg(ea);
O = CrysSym();

% ############################### Calc MisA(dg, dg_twin) ###############################
% """
% - Distance between misorientations is calculated as 1) put the two misorientations in FZ.
% 2) calc misorientation angle between two misorientations. 
% - misA formula: Kryz' thesis, P21.
% - Similar calculation: dgInFZ.m
% """

for i = 1:size(fl_obj, 1)
    g_1 = EAtoG(ea(fl_obj(i, 1), :));
    g_2 = EAtoG(ea(fl_obj(i, 2), :));

    % ----- First put the misorientation in FZ -----
    [~, ~, dg_1] = dgInFZ(g_1, g_2, O);
    
    % ----- Then treat misorientations as orientations to calculate misorientations angle (of misorientation) -----
    for idx_specialGB = 1:size(dg_obj, 3)
        dg_2 = dg_obj(:, :, idx_specialGB);
        for j = 1:24
            for k = 1:24
                gdg_1 = O(:,:,j)*dg_1;
                gdg_2 = O(:,:,k)*dg_2;

                ddg = gdg_1*gdg_2';
                ang = acosd(0.5*(trace(ddg)-1));
                if 1 - abs(0.5*(trace(ddg)-1)) < 0.001
                    ang = round(ang);
                end
                if ang < dists(i, idx_specialGB)
                    dists(i, idx_specialGB) = ang;
                end

                ddg = gdg_2*gdg_1';
                ang = acosd(0.5*(trace(ddg)-1));
                if 1 - abs(0.5*(trace(ddg)-1)) < 0.001
                    ang = round(ang);
                end
                if ang < dists(i, idx_specialGB)
                    dists(i, idx_specialGB) = ang;
                end
            end
        end
    end
%     disp([mis_ang_obj(i), dist_twin(i)])
%     disp(i);
end

end

% save('dists_full.mat', 'dists_an4_full', 'dists_an5_full');



% %%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % """
% % Check the calculation with misorientation angle of the boundary itself. 
% % """
% % ##### Load data #####
% num_neigh =  h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors').';
% neigh_list =  h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList');
% mis_ang = h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/MisorientationList');
% num_neigh(1) = [];
% 
% 
% % ##### Prepare index #####
% fl_full = zeros(size(neigh_list));
% for i = 1:size(num_neigh)
%     idx_s = sum(num_neigh(1 : i-1)) + 1;
%     idx_e = sum(num_neigh(1 : i));
%     fl_full(idx_s : idx_e) = i;
% end
% fl_full = [fl_full, neigh_list];
% 
% % ##### Get misorientation distance, as misorientation angle, to twin #####
% mask = ismember(fl_full, fl_obj, 'row');
% mis_ang_tmp = mis_ang(mask);
% idx_helper = (1:size(fl_obj, 1))';
% fl_idx_tmp = [fl_obj, idx_helper];
% fl_idx_tmp = sortrows(fl_idx_tmp, [1, 2]);
% if sum(all(fl_idx_tmp(:, 1:2) == fl_full(mask, :), 2)) ~= size(fl_obj, 1)
%     warning('calcDistFromTwin.m, index not right!')
% end
% idx_convert = sortrows([idx_helper, fl_idx_tmp(:, 3)], 2);
% mis_ang_obj = mis_ang_tmp(idx_convert(:,1));
% 
% 
% %%
% % dist_twin = dist_twin_calc;
% mask_misA60 = abs(mis_ang_obj - 60) < 5;
% mask_closetwin = abs(dists_an4_full) < 5;
% 
% disp([sum(mask_misA60), sum(mask_closetwin), sum(mask_misA60==mask_closetwin & mask_closetwin==1)]);
% 
% mask_check = (mask_closetwin & mask_misA60==0);
% tmp = [idx_helper(mask_check), mis_ang_obj(mask_check), dists_an4_full(mask_check)];








