function [is_twin, incoherence] = calcTwinsFromD3D(file, faces, eps_curv, eps_area, eps_min_ang)
% ###########################################################################
% * Input
%     - faces = [n,2] = [label_a, label_b]
%         facelabels of the faces of interest
%     - eps_
%         Filter bad triangles, this requires corresponding D3D filters to be run.
% * Output
%     - is_twin = [n,1] 
%         Indicator variable: {0, 1}
%     - incoherence = [n, 1]
%         Weighted incoherence of the grain face. NOTE only triangles with
%         good quality and do not sit on triple lines are considered. 
% * Note 
%     - This file is to be used in featurePrepOthers.m
%         A more general version is calcDistFromMisorientation.m
% ###########################################################################
% ------------------ load data for debug --------------------
% file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin3_Hsmooth.dream3d';
% eps_curv = 1;
% eps_area = 7;
% eps_min_ang = 10;
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190624_tracked_faces_full.mat', 'tracked_uniqueface_an4_full')
% faces = tracked_uniqueface_an4_full;
% clear tracked_uniqueface_an4_full
% -----------------------------------------------------------
fl = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
tri_istwin = boolean(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/TwinBoundary'))';
tri_incoherence = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/TwinBoundaryIncoherence'))';
fl = sort(fl, 2);

tri_curv =  roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
tri_area = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
tri_min_ang = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
tri_nodes = 1 + h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')';
node_types = h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType')';
tri_node_types = node_types(tri_nodes);

% """
% _full & _istwin
%     Use all triangles to check if twin, but only good quality triangles to calc incoherence.
% """
fl_full = fl(all(fl>0, 2), :);
tri_istwin = tri_istwin(all(fl>0, 2));
mask_good = all(fl>0, 2) & all(tri_node_types == 2, 2) & ...
    abs(tri_curv) < eps_curv & tri_area < eps_area & tri_min_ang > eps_min_ang;
fl_good = fl(mask_good, :);
tri_incoherence = tri_incoherence(mask_good);

is_twin = boolean(zeros(size(faces, 1), 1));
incoherence = zeros(size(faces, 1), 1);

for i = 1:length(faces)
    mask_full = (fl_full(:,1) == faces(i, 1) & fl_full(:,2) == faces(i, 2));
    if all(tri_istwin(mask_full))
        is_twin(i) = true;
    else
        if any(tri_istwin(mask_full))
            warning(['calcTwinsFromD3D.m problem, pair ', num2str(i)]);
        end
    end
    
    mask_good = (fl_good(:,1) == faces(i, 1) & fl_good(:,2) == faces(i, 2));
    incoherence(i) = sum(tri_incoherence(mask_good)) / sum(mask_good);
   
end
    
end


% ################################ Old Code ################################
% ----- prepare data -----
% fl_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
% tri_istwin_an4 = boolean(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/TwinBoundary'))';
% tri_incoherence_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/TwinBoundaryIncoherence'))';
% fl_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
% tri_istwin_an5 = boolean(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/TwinBoundary'))';
% tri_incoherence_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/TwinBoundaryIncoherence'))';
% fl_an4 = sort(fl_an4, 2);
% fl_an5 = sort(fl_an5, 2);
% 
% tri_curv_an4 =  roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
% tri_area_an4 = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
% tri_min_ang_an4 = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
% tri_curv_an5 =  roundn(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
% tri_area_an5 = roundn(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
% tri_min_ang_an5 = roundn(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
% eps_curv = 1;
% eps_area = 7;
% eps_min_ang = 10;
% mask_good_an4 = all(fl_an4>0, 2) & abs(tri_curv_an4) < eps_curv & tri_area_an4 < eps_area & tri_min_ang_an4 > eps_min_ang;
% mask_good_an5 = all(fl_an5>0, 2) & abs(tri_curv_an5) < eps_curv & tri_area_an5 < eps_area & tri_min_ang_an5 > eps_min_ang;
% 
% fl_an4 = fl_an4(mask_good_an4, :);
% tri_istwin_an4 = tri_istwin_an4(mask_good_an4);
% tri_incoherence_an4 = tri_incoherence_an4(mask_good_an4);
% fl_an5 = fl_an5(mask_good_an5, :);
% tri_istwin_an5 = tri_istwin_an5(mask_good_an5);
% tri_incoherence_an5 = tri_incoherence_an5(mask_good_an5);
% 
% % ----- get twin & incoherency -----
% istwin_an4 = boolean(zeros(size(tracked_uniqueface_an4_full, 1), 1));
% incoherence_an4 = zeros(size(tracked_uniqueface_an4_full, 1), 1);
% istwin_an5 = boolean(zeros(size(tracked_uniqueface_an5, 1), 1));
% incoherence_an5 = zeros(size(tracked_uniqueface_an5, 1), 1);
% 
% for i = 1:length(tracked_uniqueface_an4_full)
%     mask_an4 = (fl_an4(:,1) == tracked_uniqueface_an4_full(i, 1) & fl_an4(:,2) == tracked_uniqueface_an4_full(i, 2));
%     if all(tri_istwin_an4(mask_an4))
%         istwin_an4(i) = true;
%     else
%         if any(tri_istwin_an4(mask_an4))
%             warning(['D3D problem, an4, pair ', num2str(i)]);
%         end
%     end
%     incoherence_an4(i) = sum(tri_incoherence_an4(mask_an4)) / sum(mask_an4);
%  
%     mask_an5 = (fl_an5(:,1) == tracked_uniqueface_an5(i, 1) & fl_an5(:,2) == tracked_uniqueface_an5(i, 2));
%     if all(tri_istwin_an5(mask_an5))
%         istwin_an5(i) = true;
%     else
%         if any(tri_istwin_an5(mask_an5))
%             warning(['D3D problem, an5, pair ', num2str(i)]);
%         end
%     end
%     incoherence_an5(i) = sum(tri_incoherence_an5(mask_an5)) / sum(mask_an5);   
% end


