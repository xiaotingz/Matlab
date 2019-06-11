% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% - One grain in d3d can corresponds to multiple grains in simulation id_new
% due to orientation merge.  
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


% --------------------- Load data ---------------------
load('/Volumes/XIAOTING/Ni/Ni_simu_ID_An4.mat');
id_given_an4 = P;
load('/Volumes/XIAOTING/Ni/Ni_simu_ID_An5.mat');
id_given_an5 = P;
file_an4 = '/Volumes/XIAOTING/Ni/simu_An4_clean_seg.dream3d';
file_an5 = '/Volumes/XIAOTING/Ni/simu_An5_clean_seg.dream3d';
clear P

% % --------------------- 1. Find id_new ---------------------
% id_new_an4 = findNewId(id_given_an4);
% id_new_an5 = findNewId(id_given_an5);


% % --------------------- 2. Find corresp between id_given & id_new ---------------------
% corresp_simu_given_new_an4 = findIdCorrespSameStateOldNew(id_given_an4, id_new_an4);
% corresp_simu_given_new_an5 = findIdCorrespSameStateOldNew(id_given_an5, id_new_an5);


% % --------------------- 3. Find corresp between id_new_an4 and id_new_an5 ---------------------
% unique_id_given_an5 = unique(id_given_an5);
% [corresp_simu_new_45, id_newclean_an4, id_newclean_an5] = findIdCorrespSimuTwoStates(id_new_an4, id_new_an5, id_corresp_an4, id_corresp_an5, unique_id_given_an5);

% % --------------------- 4. Track id_new & id_d3d for the same state ---------------------
% corresp_d3d_simu_an4 = findCorrespSimuD3D(file_an4, id_newclean_an4);
% corresp_d3d_simu_an5 = findCorrespSimuD3D(file_an5, id_newclean_an5);


% % --------------------- 5. Track id between two D3D states ---------------------
% corresp_simu_new_45 = containers.Map(corresp_simu_new_45(:, 1), corresp_simu_new_45(:, 2));
corresp_d3d_45 = findIdCorrespD3DTwoStates(corresp_d3d_simu_an4, corresp_d3d_simu_an5, corresp_simu_new_45);


















