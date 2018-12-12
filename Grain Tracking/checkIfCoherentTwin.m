% function mask_ctwin = checkIfCoherentTwin(file, face_to_calc, thres_g, thres_n)
% ############################################################################
% * Input 
%     - face_to_calc = [n, 2], labels of the grain faces of interest.
%     - thres_norm, the threshold for angle between face_avg_norm and [111]
%  * Output
%     - mask_ctwin = [n, 1]
%  * Notes
%     - Mainly designed to check if twins are coherent or not. Use together
%     with getFaceRFvecs.m
% ############################################################################
% % ------------------ load data for debug --------------------
file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
load('181107_mig_piececorresp_comb','tracked_uniqueface_an4', 'tracked_uniqueface_an5');
load('181108_data', 'rfvecs_an4')
face_to_calc = tracked_uniqueface_an4;
rfvecs = rfvecs_an4;
thres_g = 0.06;
thres_n = 10;
clear rfvecs_an4 tracked_uniqueface_an4 tracked_uniqueface_an5
% %-----------------------------------------------------------
% ##### Load and Clean Data #####
facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_normal = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
tri_area = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'))';
EA = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles'))';
EA(1,:) = [];
EA = rad2deg(EA);
mask_inner = all(facelabel > 0, 2);
facelabel = facelabel(mask_inner, :);
tri_normal = tri_normal(mask_inner, :);
tri_area = tri_area(mask_inner, :);

% ##### Prepare Symmetry Matrix #####
O_naive = CrysSym;
n_sym = size(O_naive, 3);
O = reshape(O_naive(:), 3*3, n_sym);
O = [reshape(O(1:3, :), [], 1), reshape(O(4:6, :), [], 1), reshape(O(7:9, :), [], 1)];

% ##### Prepare Coherent Twin Normal #####
% ----- ref : https://github.com/BlueQuartzSoftware/DREAM3D/blob/develop/Source/Plugins/OrientationAnalysis/OrientationAnalysisFilters/FindGBCDMetricBased.cpp -----
fixed_normal_1 = [1,1,1]';
fixed_normal_1 = fixed_normal_1/norm(fixed_normal_1);
g_fixed = AAToG(60, [1,1,1]);
fixed_normal_2 = g_fixed' * fixed_normal_1;
area_alltwin = 0.0;
area_obj = 0.0;
ang_thres = 10.0;

% fixed_normal_1 = repmat(fixed_normal_1, 1, n_sym)';
% fixed_normal_2 = repmat(fixed_normal_2, 1, n_sym)';

% ##### First Filter Out non-twin Misorientations #####
rfvec_twin = [1,1,1]/norm([1,1,1]) * tand(60/2);
rfvec_twin = repmat(rfvec_twin, size(face_to_calc, 1), 1);
%%
% [rfvecs]  = getFaceRFvecs(file, face_to_calc);
mask_ctwin = vecnorm(rfvecs - rfvec_twin, 2, 2) < thres_g;

% ##### Then Check Plane Normals #####
% angle_diff_all_1 = [];
% angle_diff_all_2 = [];
angle_diff_list = [];
angle_diff_cnt = [];
id_list = [];

% """
% Note the symmetry matrixes are transposed this way. But the group should
% stay the same. 
% """

parfor i = 1:size(face_to_calc, 1)
    if mask_ctwin(i)
        % ----- Get triangle normals in sample frame -----
        grain_A = face_to_calc(i, 1);
        grain_B = face_to_calc(i, 2);
%         grain_A = 590;
%         grain_B = 729;
        mask_objface = (facelabel(:,1) == grain_A & facelabel(:,2) == grain_B |...
            facelabel(:,1) == grain_B & facelabel(:,2) == grain_A);
        facetri_normal = tri_normal(mask_objface, :);
        facetri_area = tri_area(mask_objface);
        angle_diffs = ones(size(facetri_normal, 1), 1)*90;

        % ----- Convert triangle normals to crystal frame and apply symmetries -----
%         g1 = reshape(G(face_to_calc(i, 1), :), 3,3)';
%         g2 = reshape(G(face_to_calc(i, 2), :), 3,3)';
        g1 = EAtoG(EA(grain_A, :));
        g2 = EAtoG(EA(grain_B, :));
    
%         n_g2_list_1 = zeros(24, 3);
%         n_g2_list_2 = zeros(24*24, 3);
        
%         % ##### normal_g2 = g2 * normal_lab ##### 
%         for j = 1:size(facetri_normal, 1)
%             normal = facetri_normal(j, :)';
%             area_alltwin = area_alltwin + facetri_area(j);
%             % ##### Min Angle ##### 
%             for k = 1:n_sym
%                 g1s = O_naive(:,:,k) * g1;
%                 normal_g1 = g1s * normal;
%                 g2s = O_naive(:,:,k) * g2;
%                 normal_g2 = g2s * normal;
% %                 n_g2_list_1(k,:) = normal_g2;
%                 % ----- dist to (111), tuple [n1, -n2] -----
%                 ang_diff_tmp1 = min(acosd(dot(normal_g1, fixed_normal_1)), acosd( - dot(normal_g1, fixed_normal_1)));
%                 ang_diff_tmp2 = min(acosd(dot(- normal_g2, fixed_normal_2)), acosd( - dot(- normal_g2, fixed_normal_2)));
%                 if ang_diff_tmp1^2 + ang_diff_tmp2^2 < ang_thres^2/2
%                     area_obj = area_obj + facetri_area(j);
% %                     angle_diffs(j) = min(ang_diff_tmp1, ang_diff_tmp2);
%                     break
%                 end
%                 if min(ang_diff_tmp1, ang_diff_tmp2) < angle_diffs(j)
%                     angle_diffs(j) = min(ang_diff_tmp1, ang_diff_tmp2);
%                 end
%                 % ----- dist to (111), tuple [-n2, n1] -----
%                 ang_diff_tmp1 = min(acosd(dot(- normal_g2, fixed_normal_1)), acosd( - dot(- normal_g2, fixed_normal_1)));
%                 ang_diff_tmp2 = min(acosd(dot(normal_g1, fixed_normal_2)), acosd(- dot(normal_g1, fixed_normal_2)));
%                 if ang_diff_tmp1^2 + ang_diff_tmp2^2 < ang_thres^2/2
%                     area_obj = area_obj + facetri_area(j);
%                     angle_diffs(j) = min(ang_diff_tmp1, ang_diff_tmp2);
%                     break
%                 end
%                 if min(ang_diff_tmp1, ang_diff_tmp2) < angle_diffs(j)
%                     angle_diffs(j) = min(ang_diff_tmp1, ang_diff_tmp2);
%                 end
%             end
%         end

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
%                     n_pair_list_2((k-1)*24+l, :) = [normal_g1', normal_g2'];
                    % ----- dist to (111), tuple [n1, -n2] -----%         
                    ang_diff_tmp1 = min(acosd(dot(normal_g1, fixed_normal_1)), acosd( - dot(normal_g1, fixed_normal_1)));
                    ang_diff_tmp2 = min(acosd(dot( - normal_g2, fixed_normal_2)), acosd( - dot( - normal_g2, fixed_normal_2)));
                    if (ang_diff_tmp1^2 + ang_diff_tmp2^2)/2 < ang_thres^2
                        area_obj = area_obj + facetri_area(j);
                    end
                    if min(ang_diff_tmp1, ang_diff_tmp2) < angle_diffs(j)
                        angle_diffs(j) = min(ang_diff_tmp1, ang_diff_tmp2)
                    end
                    % ----- dist to (111), tuple [-n2, n1] -----
                    ang_diff_tmp1 = min(acosd(dot(normal_g1, fixed_normal_1)), acosd( - dot(normal_g1, fixed_normal_1)));
                    ang_diff_tmp2 = min(acosd(dot( - normal_g2, fixed_normal_2)), acosd( - dot( - normal_g2, fixed_normal_2)));
                    if (ang_diff_tmp1^2 + ang_diff_tmp2^2)/2 < ang_thres^2
                        area_obj = area_obj + facetri_area(j);
                    end
                    if min(ang_diff_tmp1, ang_diff_tmp2) < angle_diffs(j)
                        angle_diffs(j) = min(ang_diff_tmp1, ang_diff_tmp2)
                    end
                end
            end
        end

        angle_diff_list = [angle_diff_list; angle_diffs];
        angle_diff_cnt = [angle_diff_cnt; sum(mask_objface)];
        id_list = [id_list; i];

%     % ##### check if normal_g2 = dgT * normal_g1 and normal_g2 = g2 * normal_lab give the same results ##### 
%     if sum(unique(round(n_g2_list_1, 3), 'rows') == unique(round(n_g2_list_2, 3), 'rows')) ~= 24*3
%         warning('bad normals')
%         break
%     end
        
    % ##### Check if Angle Within Threshold #####
%     if sum(angle_diffs < thres_n) < length(angle_diffs)*0.8
%         mask_ctwin(i) = false;            
%     end
    end
end

area_all = sum(tri_area);
area_obj/area_all/0.000139158
%%
figure
histogram(angle_diff_list,'Normalization','probability')
set(gca, 'FontSize', 18);
xlabel('angle between triangle normals & (111)', 'FontSize', 18)
print('all plane normal variants to (111), min angles, twins_120318','-dpng','-r300')


%% ##### Shown Coherent triangles on twin face #####
obj_facelabel_an4 = face_to_calc(i, :);
color1 = [0, 0.4470, 0.7410];
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/181107.mat','tri_node_an4', 'node_coord_an4')
tri_centr = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceCentroids'))';
tri_centr = tri_centr(mask_inner, :);

tri_connect = tri_node_an4(mask_objface, :);

figure
trisurf(tri_connect, node_coord_an4(:,1), node_coord_an4(:,2), node_coord_an4(:,3),'Facecolor',color1, 'Facealpha', 0.3, 'edgealpha', 0.3);
hold on
trisurf(tri_connect(angle_diffs < 10, :), node_coord_an4(:,1), node_coord_an4(:,2), node_coord_an4(:,3),'Facecolor',color1, 'Facealpha', 1, 'edgealpha', 1);
rotate3d on
% quiver3(tri_centr(mask_objface,1),tri_centr(mask_objface,2),tri_centr(mask_objface,3), ...
%      facetri_normal(:,1),facetri_normal(:,2),facetri_normal(:,3),1.5,'color','r');

daspect([1 1 1])

%% ##### Save Coherent Twin Figure #####
FFId = 3071;
name = ['gA_', num2str(grain_A), '_gB_', num2str(grain_B), '_FFId_', num2str(FFId), '_B']
print(name,'-dpng','-r300')



%% ##### Orientation matrix from D3D #####
file_g = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth_OrientMat.dream3d';
G = double(h5read(file_g,'/DataContainers/ImageDataContainer/CellFeatureData/AvgOrientMat'))';
% ----- G in D3D follow row first order -----
G(1, :) = [];



%% ##### Plane normal to Aenith angles #####
normal = [1,1,1];
normal = normal/norm(normal);
zenith = [atan2d(normal(2), normal(1)), 90 - acosd(normal(3))];

an4_GBCD(all(abs(an4_GBCD(:,1:2) - zenith) < 5, 2), :);

%%
% ##################################### Random Normals, Kryz #####################################
n_sample_pts = 3000;
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

figure
scatter3(sample_pts(:,1), sample_pts(:,2), sample_pts(:,3));
daspect([1 1 1])
print('uniform sampled_pts','-dpng','-r300')

% ##### Normal Distance as Angle #####
normal_ctwin_vari = reshape(O * fixed_normal_1(1,:)', 3, n_sym)';
normal_ctwin_vari = unique(normal_ctwin_vari, 'rows');
normal_ctwin_mat = normal_ctwin_vari(:);
normal_ctwin_mat = repmat(normal_ctwin_mat, 1, n_sym)';
normal_ctwin_mat = reshape(normal_ctwin_mat, [], 3);

% ang_diff_nosym = acosd(dot(sample_pts, normal_ctwin_mat, 2));
% 
% figure
% histogram(ang_diff_nosym,'Normalization','probability')
% set(gca, 'FontSize', 18);
% xlabel('uniform normals, no sym, angle distance to \{111\}', 'FontSize', 18)

ang_diffs_sym = zeros(size(sample_pts, 1), 1);
for i = 1:size(sample_pts, 1)
    normal = sample_pts(i,:)';
    normal = reshape(O * normal, 3, n_sym)';
    normal = repmat(normal, size(normal_ctwin_vari, 1), 1);
    ang_diff_tmp = acosd(dot(normal, normal_ctwin_mat, 2));
    ang_diffs_sym(i) = min(ang_diff_tmp);
end

figure
histogram(ang_diffs_sym,'Normalization','probability')
set(gca, 'FontSize', 18);
xlabel('uniform normals, sym\_min, angle distance to \{111\}', 'FontSize', 18)
print('uniform sampled_pts angle distance to {111} v2','-dpng','-r300')

%%
% ##################################### Check Mutiplication effect #####################################
n1 = [5,2,3]';
n1 = n1/norm(n1);
g1 = eye(3);
g2 = AAToG(60, [1,1,1]);

O = CrysSym;
n_sym = size(O, 3);
O = reshape(O(:), 3*3, n_sym);
O = [reshape(O(1:3, :), [], 1), reshape(O(4:6, :), [], 1), reshape(O(7:9, :), [], 1)];

gg1 = O * g1;
gg2 = O * g2;

n_vari = [];
for i = 1:24
    dg = gg1 * gg2((i-1)*3+1 : i*3, :)';
    n_vari_tmp1 = reshape(dg*n1, 3, 24)';
    dg = (gg2((i-1)*3+1 : i*3, :) * gg1')';
    n_vari_tmp2 = reshape(dg*n1, 3, 24)';
    n_vari = [n_vari; n_vari_tmp1; n_vari_tmp2];
end

unique(sort(abs(n_vari), 2), 'rows');
