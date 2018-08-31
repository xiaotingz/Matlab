% %% get delete_ID(index for grains to remove) --- remove a outer grain if they small number of faces
% % important data:
%     % data_grain [grainId, grainDiameter, #Faces, grainCurvature]
%     
% % get the Ids of Grains contacting the outer surface: by triangle label
% tmp1 = facelabel(1,:);
% tmp2 = facelabel(2,:);
% 
% boolean1 = (tmp1 == -3 | tmp1 == 0 | tmp2 == -3 | tmp2 == 0);
% data(1,:) = double(tmp1(boolean1));
% data(2,:) = double(tmp2(boolean1));
% 
% outer_raw = unique(data);
% 
% % prepare grain curvature calculation: get the Id list for inside grains
%     % delete the -3 in the outer_raw, then add 1 to every GrainID to make it index in matlab
% outer_raw(1) = [];
% tmp3 = ones(length(outer_raw),1);
% outer_ind = outer_raw + tmp3;
% 
% % if the grain contacts the outer side and has less than the threshold number of faces, delete
% delete_bool = ones(length(data),1);
% threshold = 30;
% for i = 1: length(outer_ind)
%     if num_of_neigh(outer_ind(i)) < threshold;
%         delete_bool(i) = 0;
%     end
% end
% delete_bool = logical(delete_bool);

%% get delete_ID(index for grains to remove) --- remove a outer grain if its centroid is falls within 2<R> from outer frame  
aveD_raw = sum(grain_diameter_raw) / length(grain_diameter_raw);

delete_bool = ones(length(centroids),1);
for i = 1:length(centroids)
    if centroids(i,1) < aveD_raw || X - centroids(i,1) < aveD_raw
        delete_bool(i) = 0;
    elseif centroids(i,2) < aveD_raw || Y - centroids(i,2) < aveD_raw
        delete_bool(i) = 0;
    elseif centroids(i,3) < aveD_raw || Z - centroids(i,3) < aveD_raw
        delete_bool(i) = 0;
    end
end
delete_bool = logical(delete_bool);

%% get the indexes of grains for calculation, calculate grain curvature for them
ID_list = (0:length(grain_diameter_raw)-1).';
ID_ForCal = ID_list(delete_bool);
D_ForCal = grain_diameter_raw(delete_bool);
NNeigh_ForCal = num_of_neigh(delete_bool);
% % V6
% grain_ForCal = [ID_ForCal,D_ForCal.',NNeigh_ForCal.'];
% V4
grain_ForCal = [ID_ForCal,D_ForCal,NNeigh_ForCal];

% calculate grain curvature
s_curv = 0;
s_area = 0;
tmp4 = zeros(length(grain_ForCal),4);

for i = 1:length(tmp4)
    for j = 1:length(data_final)
        if data_final(j,1) == grain_ForCal(i,1)
            s_area = s_area + data_final(j,3);
            s_curv = s_curv - data_final(j,3)*data_final(j,4);
        elseif data_final(j,2) == grain_ForCal(i,1)
            s_area = s_area + data_final(j,3);
            s_curv = s_curv + data_final(j,3)*data_final(j,4);
        end
    end
    
    tmp4(i,1) = grain_ForCal(i,1);
    tmp4(i,2) = s_area;
    tmp4(i,3) = s_curv;
    tmp4(i,4) = s_curv/s_area;
    
    s_curv = 0;
    s_area = 0;
end

% final grain data: [grainId, grainDiameter, #Faces, grainCurvature]
data_grain = [grain_ForCal,tmp4(:,3)];
