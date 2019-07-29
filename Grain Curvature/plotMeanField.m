function plotMeanField(file, data_grain, start, step, width, mean_field_id)
% ##########################################################################
% * Inputs
%     - file
%         Mainly used to read neighbor_list.
%     - data_grain = [grainId, grainDiameter, #Faces, #edges, IntegralGrainCurvature]
%         Returned by calcGrainCurvature.m 
%     - start, step, width
%         For plotting binned data. 
%     - mean_field_id
%         column_id of the mean field of interest. If not specified, use F.
% * Notes
%     - This function is to compare F-<Fnn> in G_F_mF.m USE ONLY IF no grains has
%     been filtered.
%     - This function only works for ONE VOLUME case. If two volumes, needs
%     to calculate curv_FmF separately, concatinate together and then move
%     on. 
% ##########################################################################
% ----------------------- load debug data -----------------------
% start = -80;
% width = 1;
% step = 160;
% ---------------------------------------------------------------
% ----- if mean_field not specified, use default (F) -----
if nargin == 5
    mean_field_id = 3;
end
    
ID_ForCal = data_grain(:,1);
mean_field = data_grain(:, mean_field_id);
neigh_list = h5read(file, '/DataContainers/ImageDataContainer/CellFeatureData/NeighborList');
num_faces = h5read(file, '/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors');
num_faces(1) = [];

F_mF_diff = zeros(size(ID_ForCal));
for i = 1 : length(ID_ForCal)
    grainID = ID_ForCal(i);
    total = 0;
    nlist_start = sum(num_faces(1:(grainID-1))) + 1;
    nlist_end = sum(num_faces(1:grainID));
        for j = nlist_start : nlist_end
            k = neigh_list(j);
            total = total + mean_field(k);
        end
    F_mF_diff(i) = mean_field(grainID) - total/num_faces(grainID);
end

curv_FmF = [F_mF_diff, data_grain(:,5)./(data_grain(:,2).*((4/3*pi)^(1/3)/2))];

% ------------------------------ Two volume case ------------------------------
% curv_FmF = [curv_FmF_1; curv_FmF_2];
% -----------------------------------------------------------------------------

scatter(curv_FmF(:,1),curv_FmF(:,2),10,[0.5,0.5,0.5],'filled');
hold on

data_grid = zeros(step,3);
cnt = 0;
total = 0;
for i = 1:step
% mean Curv for each Bin
    for j = 1:length(curv_FmF)
        if curv_FmF(j,1) >= start && curv_FmF(j,1) < start+width
            total = total + curv_FmF(j,2);
            cnt = cnt + 1;
        end
    end
    data_grid(i,1) = start + width/2;
    data_grid(i,2) = total/cnt;
    data_grid(i,3) = cnt;
    
    start = start + width;
    total = 0;
    cnt = 0;
end
data_grid(data_grid(:,3) <= 3,:) = [];
scatter(data_grid(:,1), data_grid(:,2), 60, 'r', 'filled', 's')

end