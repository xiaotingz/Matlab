function neighbors = getNeighList(grainID, numNeigh, NeighborList)
% ############################################################################
% NOTES
% - assuming the non-function 0 at the first element of numNeighbors has been deleted.
% ############################################################################

neighbors =  NeighborList(sum(numNeigh(1:grainID-1)) + 1 : sum(numNeigh(1:grainID)));
end