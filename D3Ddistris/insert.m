function outCell = insert(inCell, toAddList)
% ############################################################################
% NOTES
% - This file is a sub-function for findQuadPoints
% - inCell = [1, n]
% - toAddList = [1, n]
% ############################################################################
% ----------------------- load debug data -----------------------
% a = {{1,2,3};{4};{}};
% inCell = a{2};
% toAddList = [1,2];
% ---------------------------------------------------------------
inList = cell2mat(inCell);
outList = unique([inList, toAddList]);
outCell = num2cell(outList);
end

