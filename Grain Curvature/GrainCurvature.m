% important data:
    % data_grain [grainId, grainDiameter, #Faces, grainCurvature]

% get the Ids of Grains contacting the outer surface: by triangle label
tmp1 = facelabel(1,:);
tmp2 = facelabel(2,:);

boolean1 = (tmp1 == -3 | tmp1 == 0 | tmp2 == -3 | tmp2 == 0);
data(1,:) = double(tmp1(boolean1));
data(2,:) = double(tmp2(boolean1));

outer_raw = unique(data);

% prepare grain curvature calculation: get the Id list for inside grains
outer_raw(1) = [];
tmp3 = ones(length(outer_raw),1);
outer_ind = outer_raw + tmp3;
grain_diameter_raw = roundn(h5read('/Volumes/RESEARCH/Jul.7 SmallIN100/Jul7_neighborlist_qm.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
grain_Xaxis = [(0:length(grain_diameter_raw)-1).', grain_diameter_raw, num_of_neigh];
grain_Xaxis(outer_ind,:) = [];

% calculate grain curvature
s_curv = 0;
s_area = 0;
tmp4 = zeros(length(grain_Xaxis),4);

for i = 1:length(tmp4)
    for j = 1:length(data_final)
        if data_final(j,1) == grain_Xaxis(i,1)
            s_area = s_area + data_final(j,3);
            s_curv = s_curv - data_final(j,3)*data_final(j,4);
        elseif data_final(j,2) == grain_Xaxis(i,1)
            s_area = s_area + data_final(j,3);
            s_curv = s_curv + data_final(j,3)*data_final(j,4);
        end
    end
    
    tmp4(i,1) = grain_Xaxis(i,1);
    tmp4(i,2) = s_area;
    tmp4(i,3) = s_curv;
    tmp4(i,4) = s_curv/s_area;
    
    s_curv = 0;
    s_area = 0;
end

% final grain data: [grainId, grainDiameter, #Faces, grainCurvature]
data_grain = [grain_Xaxis,tmp4(:,4)];
