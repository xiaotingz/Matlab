function [x_to_y, y_to_x] = solveNodeCorresp(X, Y)
% ##########################################################################
% * Input
%     - X = [m, 3], node coordinates on the before state face
%     - Y = [n, 3], node coordinates on the after state face
% * Output
%     - x_to_y = [m, 1], for each node in X, id of the corresponding node in Y
%     - y_to_x = [n, 1], for each node in Y, id of the corresponding node in X
% * NOTE
%     - Make sure samples are aligned by changing ORIGIN.
%     - This script implements Sid's algorithm.
% ##########################################################################
% ----------------------- load debug data -----------------------
% X = center_an4;
% Y = center_an5;
% ---------------------------------------------------------------
m = size(X,1);
n = size(Y,1);

diff_tmp1 = repmat(X, 1, 1, n);
diff_tmp1 = permute(diff_tmp1, [1,3,2]);
diff_tmp2 = repmat(Y, 1, 1, m);
diff_tmp2 = permute(diff_tmp2, [3,1,2]);
diff = diff_tmp1 - diff_tmp2;
norm_diff = sum(diff .* diff, 3);

coeff = optimvar('coeff', m, n, 'LowerBound', 0);

prob = optimproblem('Objective', sum(sum(coeff .* norm_diff)));

prob.Constraints.cons1 = sum(coeff,2) == ones(m,1);
prob.Constraints.cons2 = sum(coeff,1) == (ones(1,n) * m/n);

[linsol, fval] = solve(prob);
node_corresp = linsol.coeff;


[~, x_to_y] = max(node_corresp, [], 2);
[~, y_to_x] = max(node_corresp, [], 1);

end


%%
% node_corresp_new = node_corresp;
% tmp2 = node_corresp_new(:,2);
% tmp3 = node_corresp_new(:,3);
% node_corresp_new(:,2) = node_corresp_new(:,8);
% node_corresp_new(:,3) = tmp2;
% node_corresp_new(:,8) = tmp3;
% [~, x_to_y_new] = max(node_corresp_new, [], 1);
% figure
% visualizeFace(obj_face, center_an4, center_an5, x_to_y_new)
% 
% sum(sum(node_corresp .* normDiff))
% sum(sum(corresp_new .* normDiff))