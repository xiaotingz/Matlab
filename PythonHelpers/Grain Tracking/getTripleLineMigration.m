% function [triple_line_tracked, migration] = getTripleLineMigration(tracked_tripleline, tracked_uniqueface_1, tracked_uniqueface_2, X_to_Y)
% % ############################################################################
% * Input
%     - triple_line_ = [n, 3], calculated by trackTripleLine.m
%     - tracked_uniqueface_ = [m, 2], calculated by trackUniqueFace.m
%     - X_to_Y | Y_to_X = {m, 1}, calculated in main_TrackedNodes
% * Output
%     - triple_line_tracked = [n', 3], ids for grains encapsolating the tracked triple lines 
%     - migration = [n', 1], average migration distance of the tracked triple lines
% * Notes
%     - Inputs MUST COME IN PAIR
%         (triple_line_1 & X_to_Y) || (triple_line_2 & Y_to_X) 
% % ############################################################################
% ----------------------- load debug data -----------------------
tracked_tripleline = trackTripleLine(triple_line_1, triple_line_2, look_up_table);
tracked_uniqueface_1 = tracked_uniqueface_an4;
tracked_uniqueface_2 = tracked_uniqueface_an5;
% ---------------------------------------------------------------

% ##### Get node_id for All Nodes on tracked_triple_line #####




% ##### Solve Correspondence for Migration Distance #####


