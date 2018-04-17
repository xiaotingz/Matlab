curvature_of_triangle = h5read('E:\Apr. 7\Mar.24 curvature.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures');

[r1 c1] = find(curvature_of_triangle > 10);
[r2 c2] = find(curvature_of_triangle > 50);
[r3 c3] = find(curvature_of_triangle > 100);

% [r4 c4] = find(curvature_of_triangle > 500);
% [r5 c5] = find(curvature_of_triangle > 1000);
% 
% data_forcheck = sortrows(data_final,3);
% r_max = length(data_forcheck);
% 
% id_list = [0:size(num_of_neigh)-1].';
% grainface_list = [id_list, num_of_neigh];
% grainface_sorted = sortrows(grainface_list,2);
% 
% checks.nf_total = length(x);
% checks.ntri_over100 = length(r4);
% checks.ntri_over500 = length(r2);
% checks.ntri_over1000 = length(r3);
% checks.cur_over10 = length(r4);
% checks.cur_over50 = length(r5);
% checks.maxface1 = [data_forcheck(r_max,1),data_forcheck(r_max,2)];
% checks.maxface2 = [data_forcheck(r_max-1,1),data_forcheck(r_max-1,2)];
% checks.maxface3 = [data_forcheck(r_max-2,1),data_forcheck(r_max-2,2)];
% checks.mostneighs1 = grainface_sorted((lengthgrainface_sorted)-1,1);
% checks.mostneighs2 = grainface_sorted((lengthgrainface_sorted)-2,1);
% checks.mostneighs3 = grainface_sorted((lengthgrainface_sorted)-3,1);
