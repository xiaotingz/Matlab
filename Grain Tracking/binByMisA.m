function data_grid = binByMisA(property, misAs, gridSize)
% ###########################################################################
% * misA UNIT: DEGREE!!
% * Output
%     - data_grid = [leftBoundaryOfCurrentBin, avg(property), countInBin]
% * Note
%     - The property and misAs are in the same order. Namely, an implicit correspondence exists between the two lists' element with the same index.
%     - This function bins the property by its corresponding misA and returns a histogram containing the avg(property) for each misA
% ###########################################################################
% ------------------ load data for debug --------------------
% property = face_curvDiff;
% misAs = misAs_An4(faceCorresp(:,1));
% gridSize = 0.2;
% -----------------------------------------------------------
    %   ### the max possible misA of cubic symmetry is 62.8 deg ###
    % steps = ceil(63 / stepSize);
    xAxis = (0:gridSize:63).';
    data_grid = [xAxis, zeros(length(xAxis),2)];

    %   ### the index of misA is in its value implicitly ###
    data_sorted = sortrows([misAs, property],1);
    data_sorted(:,1) = ceil(data_sorted(:,1)/gridSize);
    for i = 1 : length(data_sorted)
        ind = data_sorted(i,1);
        data_grid(ind,2) = data_grid(ind,2) + data_sorted(i,2);
        data_grid(ind,3) = data_grid(ind,3) + 1;
    end
    data_grid(:,2) = data_grid(:,2)./data_grid(:,3);
end



    



