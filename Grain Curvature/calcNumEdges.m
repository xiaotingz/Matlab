function numEdges = calcNumEdges(numNeigh, neighborList)
% ############################################################################
% NOTES
% - assuming the non-function 0 at the first element of numNeighbors has been deleted.
% ############################################################################
% numNeigh(1) = [];

numEdges = zeros(length(numNeigh), 1);
for i = 1:length(numNeigh)
    keyGrain = i;
    keyNeighbors = getNeighList(keyGrain, numNeigh, neighborList);

%     commonNeighPairs = [];
    edges = [];
    for j = 1:numNeigh(keyGrain)
        toLook = keyNeighbors(j);
        list2 = getNeighList(toLook, numNeigh, neighborList);
        intersection = Intersect(keyNeighbors, list2);
%         tmp = [int32(ones(size(intersection))).*toLook, intersection];
%         commonNeighPairs = vertcat(commonNeighPairs, tmp);
        edges = vertcat(edges, intersection);
    end

% -- Each edge is counted twice, in theory, the real number of edges should be
% exactly half of the length of our list.
    try 
        numEdges(i) = length(edges)/2; 
    catch rem(length(edges),2) ~= 0
        warning([num2str(3),' th grain, #edges ~= 2n']);
    end
end
    
end
