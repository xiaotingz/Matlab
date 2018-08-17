% function sampled_node = sampleNodeMaxMinD(node)
% ##########################################################################
% * Input
%     - node = [n, 3], the nodes on the surface
% * Output
%     - sampled_node = [m, 1], the sampled nodes on the surface
% * NOTE
%     - make sure samples are aligned by changing ORIGIN.
%     - the ideal is to make the points distribute as uniform as possible.
%     The problem is the p-dispersion problem. The objective function is
%     written as max(min(D)). See: Discrete Network Location Models by John Current
% ##########################################################################
% ----------------------- load debug data -----------------------
load('node_coord.mat');
load('4face_node_infos.mat')
obj_face = face_node_info_4;
node = node_coordsmooth_an4(obj_face{2,1},:);
num_sample = 10;
% [X, Y] = meshgrid(1:5,1:5);
% Z = zeros(5);
% node = [reshape(X,[],1), reshape(Y,[],1), reshape(Z,[],1)];
% num_sample = 10;
% ---------------------------------------------------------------
n = size(node, 1);

tmp1 = repmat(node, 1, 1, n);
tmp1 = permute(tmp1, [1,3,2]);
tmp2 = permute(tmp1, [2,1,3]);
dist = tmp1 - tmp2;
dist = vecnorm(dist, 2, 3);
% --- only the upper matrix is needed ---
dist = triu(dist);
% --- M is just a large constant ---
M = 1e5;
% --- make_matrix is to help write minDist into matrix ---
make_matrix = ones(n);
make_matrix = triu(make_matrix);
make_matrix = make_matrix - eye(n);
M = M * make_matrix;


% ##### Formulation as A Usual Optimization Problem #####
x = optimvar('x', n, 'Type','integer', 'LowerBound', 0, 'UpperBound', 1);
miist = optimvar('minDist', 'LowerBound', 0);

prob = optimproblem('Objective', - min_dist);

prob.Constraints.cons1 = sum(x) == num_sample;
prob.Constraints.cons2 = min_dist .* make_matrix + (M - dist) .* (repmat(x, 1, n)' + repmat(x, 1, n)) <= 2*M - dist;

[linsol, fval] = solve(prob);
sampled_node = linsol.x;
linsol.minDist;


% ##### Use intlinprog #####
f = [-1, zeros(n,1)];
intcon = 2:n+1;
Aeq = [0, ones(n,1)];
beq = num_sample;




% scatter(node(:,1), node(:,2), 'filled')
% hold on
% mask_centers = logical(round(sampled_node));
% centers = node(mask_centers, :);
% scatter(centers(:,1), centers(:,2), 'filled')
% end
