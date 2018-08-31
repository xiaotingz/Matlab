% ###############################################################
% criterion = 'centroidPos' | 'touchingFS' | 'numFaces' | 'NN_centoridPos' | 'NN_touchingFS'
% ###############################################################
function grain_ForCal = filterGrains(criterion, facelabel, num_of_neigh, neighborList, X,Y,Z, centroids, grain_diameter_raw)

delete_bool = ones(length(centroids),1);
try 
    if strcmp(criterion,'numFaces')   
        % get the Ids of Grains contacting the outer surface: by triangle label
        boolean1 = (facelabel(1,:) <= 0 | facelabel(2,:) <= 0);
        facelabel_freeSurf = facelabel(:,boolean1);
        outer_ind = unique(facelabel_freeSurf);
        outer_ind(outer_ind <= 0) = [];
        % if the grain contacts the outer side and has less than the threshold number of faces, delete
        threshold = 30;
        for i = 1: length(outer_ind)
            if num_of_neigh(outer_ind(i)) < threshold;
                delete_bool(i) = 0;
            end
        end
    elseif strcmp(criterion,'touchingFS') 
        boolean1 = (facelabel(1,:) <= 0 | facelabel(2,:) <= 0);
        facelabel_freeSurf = facelabel(:,boolean1);
        outer_ind = unique(facelabel_freeSurf);
        outer_ind(outer_ind <= 0) = [];

        % if the grain contacts the outer side and has less than the threshold number of faces, delete
        delete_bool(outer_ind) = 0;
    %  -- remove a outer grain if its centroid is falls within 2<R> from outer frame
    elseif strcmp(criterion,'centroidPos')
        aveD_raw = sum(grain_diameter_raw) / length(grain_diameter_raw);
        for i = 1:length(centroids)
            if centroids(i,1) < aveD_raw || X - centroids(i,1) < aveD_raw
                delete_bool(i) = 0;
            elseif centroids(i,2) < aveD_raw || Y - centroids(i,2) < aveD_raw
                delete_bool(i) = 0;
            elseif centroids(i,3) < aveD_raw || Z - centroids(i,3) < aveD_raw
                delete_bool(i) = 0;
            end
        end
    %  -- remove a outer grain if its centroid, and all its neighbors centroids, are within 2<R> from outer frame
    elseif strcmp(criterion,'NN_centoridPos')
        aveD_raw = sum(grain_diameter_raw) / length(grain_diameter_raw);
%         total = 0;
        for i = 1:length(centroids)
        % check centroid of the grain i 
            if centroids(i,1) < aveD_raw || X - centroids(i,1) < aveD_raw
                delete_bool(i) = 0;
            elseif centroids(i,2) < aveD_raw || Y - centroids(i,2) < aveD_raw
                delete_bool(i) = 0;
            elseif centroids(i,3) < aveD_raw || Z - centroids(i,3) < aveD_raw
                delete_bool(i) = 0;
            end
        % check centroid of neighbors of the ith grain
            % get neighbor list of the grain i
            nlist_start = sum(num_of_neigh(1:(i-1))) + 1;
            nlist_end = sum(num_of_neigh(1:i));
            % for the jth neighbor grain k, check its centroid
            for j = nlist_start : nlist_end
                k = neighborList(j);
%                 total = total + k;
                if centroids(k,1) < aveD_raw || X - centroids(k,1) < aveD_raw
                    delete_bool(i) = 0;
                elseif centroids(k,2) < aveD_raw || Y - centroids(k,2) < aveD_raw
                    delete_bool(i) = 0;
                elseif centroids(k,3) < aveD_raw || Z - centroids(k,3) < aveD_raw
                    delete_bool(i) = 0;
                end
            end
        end
    elseif strcmp(criterion,'NN_touchingFS') 
        boolean1 = (facelabel(1,:) <= 0 | facelabel(2,:) <= 0);
        facelabel_freeSurf = facelabel(:,boolean1);
        outer_ind = unique(facelabel_freeSurf);
        outer_ind(outer_ind <= 0) = [];

        % if the grain contacts the outer side and has less than the threshold number of faces, delete
        delete_bool(outer_ind) = 0;
        fullIDs = [1:length(delete_bool)];
        candidates = fullIDs(logical(delete_bool));
        surfGrains = fullIDs(~logical(delete_bool));
        for i = 1:length(candidates)
            neighbors = getNeighList(candidates(i), num_of_neigh, neighborList);
            if ~isempty(intersect(neighbors, surfGrains))
                delete_bool(candidates(i)) = 0;
            end
        end
    end
    if (~(strcmp(criterion,'numFaces') || strcmp(criterion,'centroidPos') || strcmp(criterion,'touchingFS') || ...
        strcmp(criterion,'NN_centoridPos') || strcmp(criterion,'NN_centoridPos') || strcmp(criterion,'NN_touchingFS')))
        warning('Wrong input criterion. Using all grains!');
    end
catch
    warning('There is a problem with the function code. Using all grains!');
end

delete_bool = logical(delete_bool);

ID_list = (1:length(grain_diameter_raw)).';
ID_ForCal = ID_list(delete_bool);
D_ForCal = grain_diameter_raw(delete_bool);
NNeigh_ForCal = num_of_neigh(delete_bool);
numEdges = calcNumEdges(num_of_neigh, neighborList);
numEdges_ForCal = numEdges(delete_bool);

% % V6
% grain_ForCal = [ID_ForCal,D_ForCal.',NNeigh_ForCal.'];
% V4
grain_ForCal = [ID_ForCal, D_ForCal, NNeigh_ForCal, numEdges_ForCal];
end