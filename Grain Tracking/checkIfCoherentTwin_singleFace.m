file = ('/Users/xiaotingzhong/Dropbox/An4new6_fixedOrigin_smooth_crop_twin_pair4524.dream3d');
grain_A = 8;
grain_B = 10;

facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_normal = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
tri_area = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'))';
% tri_centr = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceCentroids'))';
tri_node = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
EA = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles'))';
G = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/Gs'))';
EA(1,:) = [];
G(1, :) = [];
EA = rad2deg(EA);
mask_inner = all(facelabel > 0, 2);
facelabel = facelabel(mask_inner, :);
tri_normal = tri_normal(mask_inner, :);
tri_area = tri_area(mask_inner, :);
% tri_centr = tri_centr(mask_inner, :);
tri_node = tri_node(mask_inner, :);

% ##### Prepare Symmetry Matrix #####
O_naive = CrysSym;
n_sym = size(O_naive, 3);
O = reshape(O_naive(:), 3*3, n_sym);
O = [reshape(O(1:3, :), [], 1), reshape(O(4:6, :), [], 1), reshape(O(7:9, :), [], 1)];

% ##### Prepare Coherent Twin Normal #####
% ----- ref : https://github.com/BlueQuartzSoftware/DREAM3D/blob/develop/Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/FindGBCDMetricBased.cpp -----
fixed_normal_1 = [-1, 0, 0]';
fixed_normal_1 = fixed_normal_1/norm(fixed_normal_1);
dg_fixed = AAToG(60, [1,1,1]);
fixed_normal_2 = dg_fixed' * fixed_normal_1;
area_alltwin = 0.0;
area_obj = 0.0;
ang_thres = 5.0;


% ################################################################################
mask_objface = (facelabel(:,1) == grain_A & facelabel(:,2) == grain_B |...
    facelabel(:,1) == grain_B & facelabel(:,2) == grain_A);
facetri_normal = tri_normal(mask_objface, :);
facetri_area = tri_area(mask_objface);
angle_diffs = ones(size(facetri_normal, 1), 2)*90;

% ----- Convert triangle normals to crystal frame and apply symmetries -----
%         g1 = reshape(G(face_to_calc(i, 1), :), 3,3)';
%         g2 = reshape(G(face_to_calc(i, 2), :), 3,3)';
g1 = EAtoG(EA(grain_A, :));
g2 = EAtoG(EA(grain_B, :));

n_pair_list_1 = zeros(24, 6);
n_pair_list_2 = zeros(24*24, 6);

% ##### normal_g2 = g2 * normal_lab ##### 
% for j = 1:size(facetri_normal, 1)
%     normal = facetri_normal(j, :)';
%     area_alltwin = area_alltwin + facetri_area(j);
%     % ##### Min Angle ##### 
%     for k = 1:n_sym
%         g1s = O_naive(:,:,k) * g1;
%         normal_g1 = g1s * normal;
%         g2s = O_naive(:,:,k) * g2;
%         normal_g2 = g2s * normal;
%         n_pair_list_1(k,:) = [normal_g1', normal_g2'];
%         % ----- dist to (111), tuple [n1, -n2] -----
%         ang_diff_tmp1 = min(acosd(dot(normal_g1, fixed_normal_1)), acosd( - dot(normal_g1, fixed_normal_1)));
%         ang_diff_tmp2 = min(acosd(dot( - normal_g2, fixed_normal_2)), acosd( - dot( - normal_g2, fixed_normal_2)));
%         if (ang_diff_tmp1^2 + ang_diff_tmp2^2)/2 < ang_thres^2
%             area_obj = area_obj + facetri_area(j);
%         end
%         if ang_diff_tmp1 + ang_diff_tmp2 < sum(angle_diffs(j, :))
%             angle_diffs(j, :) = [ang_diff_tmp1, ang_diff_tmp2];
%         end
%         % ----- dist to (111), tuple [-n2, n1] -----
%         ang_diff_tmp1 = min(acosd(dot( - normal_g2, fixed_normal_1)), acosd( - dot( - normal_g2, fixed_normal_1)));
%         ang_diff_tmp2 = min(acosd(dot(normal_g1, fixed_normal_2)), acosd( - dot(normal_g1, fixed_normal_2)));
%         if (ang_diff_tmp1^2 + ang_diff_tmp2^2)/2 < ang_thres^2
%             area_obj = area_obj + facetri_area(j);
%         end
%         if ang_diff_tmp1 + ang_diff_tmp2 < sum(angle_diffs(j, :))
%             angle_diffs(j, :) = [ang_diff_tmp1, ang_diff_tmp2];
%         end
%     end
% end

% ##### normal_g2 = dgT * normal_g1 ##### 
for j = 1:size(facetri_normal, 1)
    normal = facetri_normal(j, :)';
    area_alltwin = area_alltwin + facetri_area(j);
    
    % ##### Min Angle ##### 
    for k = 1:n_sym
        g1s = O_naive(:,:,k) * g1;
        normal_g1 = g1s * normal;
        % ----- dist to (111) -----
%                 ang_diff_tmp1 = acosd(dot(normal_g1, fixed_normal_1, 2));
%                 if min(ang_diff_tmp1) < angle_diffs(j)
%                     angle_diffs(j) = min(ang_diff_tmp1);
%                 end
        for l = 1:n_sym
            dgs = g1s*g2'*O_naive(:,:,l)';
            normal_g2 = dgs' * normal_g1;
            n_pair_list_2((k-1)*24+l, :) = [normal_g1', normal_g2'];
            % ----- dist to (111), tuple [n1, -n2] -----%         
            ang_diff_tmp1 = min(acosd(dot(normal_g1, fixed_normal_1)), acosd( - dot(normal_g1, fixed_normal_1)));
            ang_diff_tmp2 = min(acosd(dot( - normal_g2, fixed_normal_2)), acosd( - dot( - normal_g2, fixed_normal_2)));
            if (ang_diff_tmp1^2 + ang_diff_tmp2^2)/2 < ang_thres^2
                area_obj = area_obj + facetri_area(j);
            end
            % ----- dist to (111), tuple [-n2, n1] -----
            ang_diff_tmp1 = min(acosd(dot(normal_g1, fixed_normal_1)), acosd( - dot(normal_g1, fixed_normal_1)));
            ang_diff_tmp2 = min(acosd(dot( - normal_g2, fixed_normal_2)), acosd( - dot( - normal_g2, fixed_normal_2)));
            if (ang_diff_tmp1^2 + ang_diff_tmp2^2)/2 < ang_thres^2
                area_obj = area_obj + facetri_area(j);
            end
        end
    end

end


area_total = sum(tri_area);
disp(['fixedNormal = ', num2str(round(fixed_normal_1(1), 2)), ' ', num2str(round(fixed_normal_1(2), 2)), ' ', ...
    num2str(round(fixed_normal_1(3), 2)), ';  MRD = ', num2str(area_obj/area_total/0.000139158)]);


color1 = [0, 0.4470, 0.7410];
tri_connect = tri_node(mask_objface, :);

% figure
% trisurf(tri_connect, node_coord(:,1), node_coord(:,2), node_coord(:,3),'Facecolor',color1, 'Facealpha', 0.3, 'edgealpha', 0.3);
% hold on
% trisurf(tri_connect(any(angle_diffs < 10), :), node_coord(:,1), node_coord(:,2), node_coord(:,3),'Facecolor',color1, 'Facealpha', 1, 'edgealpha', 1);
% rotate3d on
% daspect([1 1 1])


%% ######################################################################################################
% """
% Sample random points on a hemisphere.
% Note the 3000 sampled points are the same as the points returned by
% Kryz's code except for the added points at the end.
% """
n_sample_pts = 6000;
sample_pts = zeros(n_sample_pts, 3);
inc = pi * (3 - sqrt(5));
off = 2.0 / n_sample_pts;
for i = 1:n_sample_pts
    y = i * off - 1.0 + 0.5*off;
    r = sqrt(max(0.0, 1.0 - y^2));
    phi = i * inc;
    z = sin(phi) * r;
    
    if z > 0
        sample_pts(i, 1) = cos(phi)*r;
        sample_pts(i, 2) = y;
        sample_pts(i, 3) = z;
    end
end
sample_pts(sum(sample_pts,2)==0, :) = [];   
% ----- Zenith angles of the sampled points -----
zenith = [atan2d(sample_pts(:,2), sample_pts(:,1)), 90 - acosd(sample_pts(:,3))];

% ##### Find Index for Direction of Interest #####
obj_direction = [1, 0, 1]';
obj_direction = obj_direction/norm(obj_direction);

[~, idx] = min(sum(abs(sample_pts - repmat(obj_direction', length(sample_pts), 1)),2));
disp(['an4,  ', '[', num2str(round(obj_direction(1),2)), ',', num2str(round(obj_direction(2),2)), ','...
    , num2str(round(obj_direction(3),2)), '],  ','metric-based_GBCD = ', num2str(an4_GBCD(idx, 3))]);


%% ######################################################################################################
% file = ('/Users/xiaotingzhong/Dropbox/An4new6_fixedOrigin_smooth_crop_twin_pair4524.dream3d');
file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';

facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_normal = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
tri_area = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'))';
tri_node = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
EA = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles'))';
istwin = boolean(h5read(file, '/DataContainers/TriangleDataContainer/FaceData/TwinBoundary'))';
EA(1,:) = [];
EA = rad2deg(EA);
% mask_inner = all(facelabel > 0, 2);
% facelabel = facelabel(mask_inner, :);
% tri_normal = tri_normal(mask_inner, :);
% tri_area = tri_area(mask_inner, :);
% tri_node = tri_node(mask_inner, :);
facelabel = facelabel(istwin, :);
tri_normal = tri_normal(istwin, :);
tri_area = tri_area(istwin, :);
tri_node = tri_node(istwin, :);

% ##### Prepare Symmetry Matrix #####
O_naive = CrysSym;
n_sym = size(O_naive, 3);
O = reshape(O_naive(:), 3*3, n_sym);
O = [reshape(O(1:3, :), [], 1), reshape(O(4:6, :), [], 1), reshape(O(7:9, :), [], 1)];

% ##### Prepare Coherent Twin Normal #####
% ----- ref : https://github.com/BlueQuartzSoftware/DREAM3D/blob/develop/Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/FindGBCDMetricBased.cpp -----
% fixed_normal_1 = sample_pts(1,:)';
fixed_normal_1 = [1, 0, 1]';
fixed_normal_1 = fixed_normal_1/norm(fixed_normal_1);
dg_fixed = AAToG(60, [1,1,1]);
fixed_normal_2 = dg_fixed' * fixed_normal_1;
area_obj = 0.0;
ang_thres = 5.0;
misa_res = 5.0;
min_angles = ones(size(facelabel, 1),1) * 360;
min_angles_2 = [];

%%
% ########################################
disp(' ');
disp('initializing ');
for i = idxes(775)
    % ----- Show process -----    
    if rem(i, length(facelabel) * 0.1) == 0
        progress = floor(i/(length(facelabel) * 0.1))*10;
        clear_info = repmat('\b', 1, 13+length(num2str(progress)));
        fprintf(clear_info);
        disp(['processed ', num2str(progress), ' %']);
    end

    % ----- Convert triangle normals to crystal frame and apply symmetries -----
    g1 = EAtoG(EA(facelabel(i, 1), :));
    g2 = EAtoG(EA(facelabel(i, 2), :));
    
    
    % ##### normal_g2 = dgT * normal_g1 ##### 
    normal = tri_normal(i, :)';

    % ##### Min Angle ##### 
    for k = 1:n_sym
        g1s = O_naive(:,:,k) * g1;
        normal_g1 = g1s * normal;
        for l = 1:n_sym
            dgs = g1s*(O_naive(:,:,l)*g2)';
            diff_fixed_misa_1 = acosd((trace(dgs * dg_fixed') - 1.0) * 0.5);
            diff_fixed_misa_2 = acosd((trace(dgs' * dg_fixed') - 1.0) * 0.5);
            if diff_fixed_misa_1 < misa_res
                normal_g2 = dgs' * normal_g1;
                for sign = 1:2
                    if sign > 1
                        sign = -1;
                    end
                    % ----- dist to (111), tuple [n1, -n2] -----%         
                    ang_diff_tmp1 = acosd(sign*(dot(normal_g1, fixed_normal_1)));
                    ang_diff_tmp2 = acosd(-sign*(dot( -normal_g2, fixed_normal_2)));
                    if (ang_diff_tmp1^2 + ang_diff_tmp2^2)/2 < ang_thres^2
                        area_obj = area_obj + tri_area(i);
                        min_angles_2 = [min_angles_2; min(ang_diff_tmp1, ang_diff_tmp2)];
                    end
%                     disp([num2str(deg2rad(ang_diff_tmp1)),',  ', num2str(deg2rad(ang_diff_tmp2))]);
%                     disp(['n_g1 = ', num2str(round(normal_g1(1),2)), ', ', num2str(round(normal_g1(2),2)), ', ', num2str(round(normal_g1(3),2))]);
%                     disp(['n_g2 = ', num2str(round(- normal_g2(1),2)), ', ', num2str(round(- normal_g2(2),2)), ', ', num2str(round(- normal_g2(3),2))]);
                    if min(ang_diff_tmp1, ang_diff_tmp2) < min_angles(i)
                        min_angles(i) = min(ang_diff_tmp1, ang_diff_tmp2);
                    end
                end
            elseif  diff_fixed_misa_2 < misa_res
                normal_g2 = dgs' * normal_g1;
                for sign = 1:2
                    if sign > 1
                        sign = -1;
                    end
                    % ----- dist to (111), tuple [n1, -n2] -----%         
                    ang_diff_tmp1 = acosd(sign*(dot( -normal_g2, fixed_normal_1)));
                    ang_diff_tmp2 = acosd(-sign*(dot(normal_g1, fixed_normal_2)));
                    if (ang_diff_tmp1^2 + ang_diff_tmp2^2)/2 < ang_thres^2
                        area_obj = area_obj + tri_area(i);
                        min_angles_2 = [min_angles_2; min(ang_diff_tmp1, ang_diff_tmp2)];
                    end
%                     disp([num2str(deg2rad(ang_diff_tmp1)),',  ', num2str(deg2rad(ang_diff_tmp2))]);
%                     disp(['n_g1 = ', num2str(round(normal_g1(1),2)), ', ', num2str(round(normal_g1(2),2)), ', ', num2str(round(normal_g1(3),2))]);
%                     disp(['n_g2 = ', num2str(round(- normal_g2(1),2)), ', ', num2str(round(- normal_g2(2),2)), ', ', num2str(round(- normal_g2(3),2))]);
                    if min(ang_diff_tmp1, ang_diff_tmp2) < min_angles(i)
                        min_angles(i) = min(ang_diff_tmp1, ang_diff_tmp2);
                    end
                end
            end
        end
    end
    
end

% area_total = sum(tri_area);
% disp(['fixedNormal = ', num2str(round(fixed_normal_1(1), 2)), ' ', num2str(round(fixed_normal_1(2), 2)), ' ', ...
%     num2str(round(fixed_normal_1(3), 2)), ';  MRD = ', num2str(area_obj/area_total/0.000139158)]);
 
%%
set(0,'defaultAxesLabelFontSize',1.2)
set(0,'defaultAxesFontSize',15)
subplot(3,1,1)
histogram(min_angles_2)
xlabel('distance to (1,0,1), added triangles')
subplot(3,1,2)
histogram(min_angles(min_angles<5))
xlabel('distance to (1,0,1), qualified triangles')
subplot(3,1,3)
histogram(min_angles(min_angles~=360))
xlabel('distance to (1,0,1), all triangles')
% print('an4_101_min_angles','-dpng','-r300')
% save('181211_an4_101', 'min_angles_2', 'min_angles' );

% idxes = 1:length(mask_inner);
% idxes = idxes(mask_inner)';
% result = [idxes(min_angles~=360), min_angles(min_angles~=360)];
% d3d_clean = d3d_dist(d3d_dist(:,2)~=360, 1:2);
% diff = result(:,2) - d3d_clean(:,2);
% sum(diff > 0.01)

% color1 = [0, 0.4470, 0.7410];
% tri_connect = tri_node(mask_objface, :);
% 
% figure
% trisurf(tri_connect, node_coord(:,1), node_coord(:,2), node_coord(:,3),'Facecolor',color1, 'Facealpha', 0.3, 'edgealpha', 0.3);
% hold on
% trisurf(tri_connect(any(angle_diffs < 10), :), node_coord(:,1), node_coord(:,2), node_coord(:,3),'Facecolor',color1, 'Facealpha', 1, 'edgealpha', 1);
% rotate3d on
% daspect([1 1 1])
