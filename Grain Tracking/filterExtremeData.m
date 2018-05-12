function filteredData = filterExtremeData(inputData, keepRatio)
% ##########################################################################
% filter the plot data by get rid of the data points of the largest Xs and Ys
% * the final keepRatio = 1 - (1 - keepRatio)*2
% ##########################################################################
    n = length(inputData);
    tmp_outliers = [(1:n)', inputData];
    tmp_outliers = abs(tmp_outliers);
    mask_outliers = ones(n, 1);
    tmp_outliers = sortrows(tmp_outliers, [2,3]);
    mask_outliers(tmp_outliers(round(keepRatio*n):n, 1)) = 0;
    tmp_outliers = sortrows(tmp_outliers, [3,2]);
    mask_outliers(tmp_outliers(round(keepRatio*n):n, 1)) = 0;
    mask_outliers = logical(mask_outliers);
    filteredData = inputData(mask_outliers,:);
end





