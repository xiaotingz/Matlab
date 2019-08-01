function dists = calcTriDistFromTargetNormals(tri_normals, tri_fls, ea, target_normals)
% ###########################################################################
% * Input
%     - tri_normals = [n,3]
%         Trianlge normal, n = #triangles
%     - tri_fls = [n, 2]
%         Facelabel
%     - ea = [k, 3]
%         Euler angles, k = #grains
%     - target_normals = [m, 3]
%         Plane normals of interest, e.g. [1,1,1], [1,1,0], [1,0,0]
% * Output
%     - dists = [n,1] 
%         Distance from a triangle normal to a target normal direction.
% * Note 
%     - This function is to be used in calcDistFromTargetNormals.m
%     - This function calculates the distance between (tri_normal).
%     - Note that distances to different target_normals can not be cross
%     compared because they have different symmetric copies. High symmetry
%     directions tend to have fewer small distances.
% ###########################################################################
% ------------------ load data for debug --------------------
% file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
% eps_curv = 1;
% eps_area = 7;
% eps_min_ang = 10;
% target_normals = [1,1,1;1,1,0;1,0,0;1,1,2];
% target_normals = target_normals ./ vecnorm(target_normals, 2, 2);
% 
% tri_normals = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/FaceNormals').';
% tri_fls = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels').';
% ea =  h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles').';
% ea(1, :) = [];
% ea = rad2deg(ea);
% 
% tri_nodes = 1 + h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')';
% node_types = h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType')';
% tri_curv = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
% tri_area = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
% tri_min_ang = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
% tri_node_types = node_types(tri_nodes);
% mask_good = all(tri_fls > 0, 2) & all(tri_node_types == 2, 2) & ... 
%     abs(tri_curv) < eps_curv & tri_area < eps_area & tri_min_ang > eps_min_ang;
% clear tri_area tri_curv tri_min_ang tri_node_types tri_nodes 
% 
% tri_fls = tri_fls(mask_good, :);
% tri_normals = tri_normals(mask_good, :);
% -----------------------------------------------------------
O = CrysSym();
dists = ones(size(tri_normals, 1), size(target_normals, 1)) * 360;

for i = 1:size(tri_normals, 1)
    n = tri_normals(i, :);
    g_1 = EAtoG(ea(tri_fls(i, 1), :));
    g_2 = EAtoG(ea(tri_fls(i, 2), :));
    
    % """
    % angle between vectors:
    %   https://www.mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab
    % ref paper for symmetries:
    %   The distribution of internal interfaces in  polycrystals
    % """
    for idx_target_normal = 1:size(target_normals, 1)
        n_obj = target_normals(idx_target_normal, :);
        for j = 1:24
            for inverse = 0:1
                gg_1 = O(:,:,j)*g_1;
                gg_2 = O(:,:,j)*g_2;

                nn_1 = gg_1*n';
                if inverse == 1
                    nn_1 = - nn_1;
                end
                ang = atan2d(norm(cross(n_obj, nn_1)),dot(n_obj, nn_1));
    %                 if 1 - abs(0.5*(trace(nn_1)-1)) < 0.001
    %                     ang = round(ang);
    %                 end
                if ang < dists(i, idx_target_normal)
                    dists(i, idx_target_normal) = ang;
                end

                nn_2 = gg_2*n';
                if inverse == 1
                    nn_2 = - nn_2;
                end
                ang = atan2d(norm(cross(n_obj, nn_2)),dot(n_obj, nn_2));
                if ang < dists(i, idx_target_normal)
                    dists(i, idx_target_normal) = ang;
                end
            end
        end
    end
end

end
%% ############################ Checks ############################
% % """
% % Idea: Compare dists to IPFColors by D3D, see if match.
% %     Small distance <-> the corresponding IPF color big and the other
% %     two IPF colors small in either crystal frame.
% % Note tri_ipf = [[100], [110], [111]]
% % """

% tri_ipf = h5read(file,'/DataContainers/TriangleDataContainer/FaceData/IPFColors')';
% tri_ipf = tri_ipf(mask_good, :);
% 
% mask_111 = dists(:,1) < 5;
% data_111 = [dists(mask_111, 1), max(tri_ipf(mask_111, 3), tri_ipf(mask_111, 6))];
% data_111_full =  [dists(mask_111, 3), dists(mask_111, 2), dists(mask_111, 1), ...
%     zeros(sum(mask_111), 1), tri_ipf(mask_111, :)];
% 
% mask_110 = dists(:,2) < 5;
% data_110 = [dists(mask_110, 2), max(tri_ipf(mask_110, 2), tri_ipf(mask_110, 5))];
% data_110_full =  [dists(mask_110, 3), dists(mask_110, 2), dists(mask_110, 1), ...
%     zeros(sum(mask_110), 1), tri_ipf(mask_110, :)];
% 
% mask_100 = dists(:,3) < 5;
% data_100 = [dists(mask_100, 3), max(tri_ipf(mask_100, 1), tri_ipf(mask_100, 4))];
% data_100_full =  [dists(mask_100, 3), dists(mask_100, 2), dists(mask_100, 1), ...
%     zeros(sum(mask_100), 1), tri_ipf(mask_100, :)];


% ----- Confirming (111) has the largest population -----
% """
% Note tri_dists can not be compared the same way because [111], [110],
%   [100] and [211] have different #(symmetry copies).
% """
% tri_ipf = int32(tri_ipf);
% 
% ipf_dists = - ones(size(tri_ipf,1), 3);
% ipf_dists(:, 1) = max(tri_ipf(:,1) - tri_ipf(:,2) - tri_ipf(:,3), ...
%                         tri_ipf(:,4) - tri_ipf(:,5) - tri_ipf(:,6));
% ipf_dists(:, 2) = max(tri_ipf(:,2) - tri_ipf(:,1) - tri_ipf(:,3), ...
%                         tri_ipf(:,5) - tri_ipf(:,4) - tri_ipf(:,6));
% ipf_dists(:, 3) = max(tri_ipf(:,3) - tri_ipf(:,2) - tri_ipf(:,1), ...
%                         tri_ipf(:,6) - tri_ipf(:,5) - tri_ipf(:,5));
% 
% sum(ipf_dists > 200)




