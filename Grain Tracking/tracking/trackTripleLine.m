function [tracked_tripleline_1, tracked_tripleline_2] = trackTripleLine(triple_line_1, triple_line_2, look_up_table)
% % ############################################################################
% * Input
%     - triple_line_ = [n, 3], calculated by findTripleLines.m in Topologies folder
%     - look_up_table, order according to an5
% * Output
%     - triple_line_tracked_ = [n', 3], ids for grains encapsolating the tracked triple lines 
% * Notes
%     - free surface contaction is not considered in this tracking
% % ############################################################################
% ----------------------- load debug data -----------------------
% triple_line_1 = triple_line_an4; 
% triple_line_2 = triple_line_an5; 
% ---------------------------------------------------------------

% ##### Make a Full Look-up Table ##### 
look_up_table = sortrows(look_up_table, 2);
look_up_2to1 = zeros(max(look_up_table(:,2)),1);
idx = 1;
for i = 1 : max(look_up_table(:,2))
    if look_up_table(idx,2) == i
        look_up_2to1(i) = look_up_table(idx,1);
        idx = idx + 1;
    else
        look_up_2to1(i) = NaN;
    end
end

look_up_table = sortrows(look_up_table, 1);
look_up_1to2 = zeros(max(look_up_table(:,1)),1);
idx = 1;
for i = 1 : max(look_up_table(:,1))
    if look_up_table(idx,1) == i
        look_up_1to2(i) = look_up_table(idx,2);
        idx = idx + 1;
    else
        look_up_1to2(i) = NaN;
    end
end

% ----- clean up and sort the list for comparison -----
triple_line_1in2 = look_up_1to2(triple_line_1);
triple_line_1in2(any(isnan(triple_line_1in2), 2), :) = [];
triple_line_1in2 = sort(triple_line_1in2, 2);
triple_line_1in2 = sortrows(triple_line_1in2);
triple_line_2 = sortrows(triple_line_2);

% ##### Tracked = Common Index #####
tracked_tripleline_2 = intersect(triple_line_2, triple_line_1in2, 'rows');
tracked_tripleline_1 = look_up_2to1(tracked_tripleline_2);

% % ##### Check #####
% tmp = look_up_1to2(triple_line_tracked_1);
% sum(triple_line_tracked_2 == tmp)

end
