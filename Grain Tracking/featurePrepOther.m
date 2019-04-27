% ##########################################################################
% * Notes
%     - Misorientaion distance is calculated from misorientation angle. 
%     - Distance to coherent twin is calculated as defined by Moraweic.
% ##########################################################################

% load('/Users/xiaotingzhong/Documents/Matlab/Grain Tracking/data/181108_rfvec.mat', 'rfvecs_an4', 'rfvecs_an5');
file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
file_an5 = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth.dream3d';
load('look_up_table_an4_an5crop.mat')

% ############### get faces to track ###############
[tracked_uniqueface_an4, tracked_uniqueface_an5] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_complete_faces');

% % ############### identify broken faces ###############
% is_onepiece = ones(size(tracked_uniqueface_an4, 1), 1);
% is_onepiece(face_piecewise) = 0;


% ############### identify twins, from D3D ###############
% ----- prepare data -----
fl_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
tri_istwin_an4 = boolean(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/TwinBoundary'))';
tri_incoherence_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/TwinBoundaryIncoherence'))';
fl_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
tri_istwin_an5 = boolean(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/TwinBoundary'))';
tri_incoherence_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/TwinBoundaryIncoherence'))';
fl_an4 = sort(fl_an4, 2);
fl_an5 = sort(fl_an5, 2);

tri_curv_an4 =  roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
tri_area_an4 = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
tri_min_ang_an4 = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
tri_curv_an5 =  roundn(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
tri_area_an5 = roundn(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
tri_min_ang_an5 = roundn(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;
mask_good_an4 = all(fl_an4>0, 2) & abs(tri_curv_an4) < eps_curv & tri_area_an4 < eps_area & tri_min_ang_an4 > eps_min_ang;
mask_good_an5 = all(fl_an5>0, 2) & abs(tri_curv_an5) < eps_curv & tri_area_an5 < eps_area & tri_min_ang_an5 > eps_min_ang;

fl_an4 = fl_an4(mask_good_an4, :);
tri_istwin_an4 = tri_istwin_an4(mask_good_an4);
tri_incoherence_an4 = tri_incoherence_an4(mask_good_an4);
fl_an5 = fl_an5(mask_good_an5, :);
tri_istwin_an5 = tri_istwin_an5(mask_good_an5);
tri_incoherence_an5 = tri_incoherence_an5(mask_good_an5);

% ----- get twin & incoherency -----
istwin_an4 = boolean(zeros(size(tracked_uniqueface_an4, 1), 1));
incoherence_an4 = zeros(size(tracked_uniqueface_an4, 1), 1);
istwin_an5 = boolean(zeros(size(tracked_uniqueface_an5, 1), 1));
incoherence_an5 = zeros(size(tracked_uniqueface_an5, 1), 1);

%%
for i = 1:length(tracked_uniqueface_an4)
    mask_an4 = (fl_an4(:,1) == tracked_uniqueface_an4(i, 1) & fl_an4(:,2) == tracked_uniqueface_an4(i, 2));
    if all(tri_istwin_an4(mask_an4))
        istwin_an4(i) = true;
    else
        if any(tri_istwin_an4(mask_an4))
            warning(['D3D problem, an4, pair ', num2str(i)]);
        end
    end
    incoherence_an4(i) = sum(tri_incoherence_an4(mask_an4)) / sum(mask_an4);
 
    mask_an5 = (fl_an5(:,1) == tracked_uniqueface_an5(i, 1) & fl_an5(:,2) == tracked_uniqueface_an5(i, 2));
    if all(tri_istwin_an5(mask_an5))
        istwin_an5(i) = true;
    else
        if any(tri_istwin_an5(mask_an5))
            warning(['D3D problem, an4, pair ', num2str(i)]);
        end
    end
    incoherence_an5(i) = sum(tri_incoherence_an5(mask_an5)) / sum(mask_an5);   
end
    

%%
% ############### distance from twins ###############
rfvec_twin = [1,1,1]/norm([1,1,1]) * tand(60/2);

% """ Note 'AvgEAs' | 'AvgEulerAngles' """
[rfvecs_an4]  = getFaceRFvecs(file_an4, tracked_uniqueface_an4);
[rfvecs_an5]  = getFaceRFvecs(file_an5, tracked_uniqueface_an5);

dist_twin_an4 = vecnorm(rfvecs_an4 - rfvec_twin, 2, 2);
dist_twin_an5 = vecnorm(rfvecs_an5 - rfvec_twin, 2, 2);


%% #################################### Write txt file ####################################
not_twin_an4 = ~istwin_an4;
not_twin_an5 = ~istwin_an5;

fileID = fopen('190425_features_otherinfo.txt','w');
% fprintf(fileID,'%s, %s, %s, %s\n','is_onepiece', 'not_twin', 'dist_twin', 'weighteddist_ctwin');
fprintf(fileID,'%s, %s, %s, %s, %s, %s\n','not_twin_an4', 'rfvect_disttwin_an4', 'incoherence_an4','not_twin_an5', 'rfvect_disttwin_an5', 'incoherence_an5');
for i = 1:length(not_twin_an4)
    fprintf(fileID, '%6d, %6.3f, %6.3f, %6d, %6.3f, %6.3f\n', ...
        not_twin_an4(i), dist_twin_an4(i), incoherence_an4(i), not_twin_an5(i), dist_twin_an5(i), incoherence_an5(i)) ;
end
fclose(fileID);

