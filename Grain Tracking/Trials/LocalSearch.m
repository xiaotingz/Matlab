% function obj_value = LocalSearch(i, nodes_obj_face, nodes_ref_face)
% ##########################################################################
% * Input
%     - i, the node whose neighborhood will be explored
%     - nodes_obj_face = [n, 3], coordinates of nodes that needs adjustment
%     - nodes_ref_face = [m, 3], coordinates of nodes on the reference face
%                   !!!!!!!!    TRY Parameter P    !!!!!!!!
% * NOTE
%     - the objective is to 1). keep the nodes in nodes_move as dispersed
%     as possible and 2). reduce the distance between nodes_move and
%     nodes_ref
% ##########################################################################
% ----------------------- load debug data -----------------------
nodes_obj_face = nodes_kmeans_an4;
nodes_ref_face = nodes_kmeans_an5;
% ---------------------------------------------------------------










