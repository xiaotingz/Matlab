% ##### Prepare Coherent Twin Normal #####
O = CrysSym;
normal_ctwin = [1,1,1]';
normal_ctwin = normal_ctwin/sqrt(3);
rfvec_twin = normal_ctwin' * tand(60/2);
normal_ctwin_vari2 = zeros(24, 3);
for i = 1:size(O, 3)
    normal_ctwin_vari2(i, :) = (O(:,:,i) * normal_ctwin)';
end
normal_ctwin_vari2 = unique(normal_ctwin_vari2, 'rows');

%%
% #######################################################################################  
% V1: normal_vari, always in crystal frame.
for j = 1:size(facetri_normal, 1)
            normal = facetri_normal(j, :)';
            % ##### Min Angle ##### 
            for l = 1:size(O, 3)
                normal_vari = (O(:,:,l) * g1) * normal;
                normal_vari = repmat(normal_vari, 1, 8)';
                ang_diff_1 = acosd(dot(normal_vari, normal_ctwin_vari, 2));
                if min(ang_diff_tmp) < angle_diffs(j)
                    angle_diffs(j) = min(ang_diff_tmp);
                end
                normal_vari = (O(:,:,l) * g2) * normal;
                normal_vari = repmat(normal_vari, 1, 8)';
                
                ang_diff_tmp = acosd(dot(normal_vari, normal_ctwin_vari, 2));
                if min(ang_diff_tmp) < angle_diffs(j)
                    angle_diffs(j) = min(ang_diff_tmp);
                end
            end
%             % ##### Kryz's way #####
%             for l = 1:size(O, 3)
%                 normal_g1 = (O(:,:,l) * g1) * normal;
%                 normal_g1 = repmat(normal_g1, 1, 8)';
% %                 ang_diff_tmp = abs(atan2d(vecnorm(cross(normal_vari, normal_ctwin_vari), 2, 2), dot(normal_vari, normal_ctwin_vari, 2)));
%                 ang_diff_1 = abs(acosd(dot(normal_g1, normal_ctwin_vari, 2)));
%                 normal_g2 = (O(:,:,l) * g2) * normal;
%                 normal_g2 = repmat(normal_g2, 1, 8)';
% %                 ang_diff_tmp = abs(atan2d(vecnorm(cross(normal_vari, normal_ctwin_vari), 2, 2), dot(normal_vari, normal_ctwin_vari, 2)));
%                 ang_diff_2 = abs(acosd(dot(normal_g2, normal_ctwin_vari, 2)));
%                 if ang_diff_1^2 + ang_diff_2^2 < thres_n^2
%                     
%                 end
%             end
        end
        
        angle_diff_list = [angle_diff_list; angle_diffs];
        angle_diff_cnt = [angle_diff_cnt; sum(mask_objface)];
        id_list = [id_list; i];
        
        % ##### Check if Angle Within Threshold #####
        if sum(angle_diffs < thres_n) < length(angle_diffs)*0.8
            mask_ctwin(i) = false;            
        end

        
% #######################################################################################  
% V2: normal_vari, n1 = n, n2 = dg'*n1, slow for loop
for j = 1:size(facetri_normal, 1)
    normal = facetri_normal(j, :)';
    % ##### Min Angle ##### 
    for l = 1:size(O, 3)
        for k = 1:size(O, 3)
            gg1 = O(:,:,l) * g1;
            gg2 = O(:,:,k) * g2;
            dg = gg1 * gg2'; 
            normal_vari = gg1 * normal;
            normal_vari = repmat(normal_vari, 1, 8)';
            ang_diff_1 = acosd(dot(normal_vari, normal_ctwin_vari, 2));
            if min(ang_diff_tmp) < angle_diffs(j)
                angle_diffs(j) = min(ang_diff_tmp);
            end
            normal_vari = gg2 * normal;
            normal_vari = repmat(normal_vari, 1, 8)';
            normal_ctwin_vari = (dg * normal_ctwin_vari')';
            ang_diff_tmp = acosd(dot(normal_vari, normal_ctwin_vari, 2));
            if min(ang_diff_tmp) < angle_diffs(j)
                angle_diffs(j) = min(ang_diff_tmp);
            end
            normal_ctwin_vari = (dg' * normal_ctwin_vari')';
            ang_diff_tmp = acosd(dot(normal_vari, normal_ctwin_vari, 2));
            if min(ang_diff_tmp) < angle_diffs(j)
                angle_diffs(j) = min(ang_diff_tmp);
            end
        end
    end
end



% #######################################################################################  
% V3: normal, n1 = g1*n, n2 = g2*n
angle_diff_list = []
for i = 1:length(face_to_calc)
    disp(i)
    if mask_ctwin(i)
        % ----- Get triangle normals in sample frame -----
        grain_A = face_to_calc(i, 1);
        grain_B = face_to_calc(i, 2);
%         grain_A = 590;
%         grain_B = 729;
        mask_objface = (facelabel(:,1) == grain_A & facelabel(:,2) == grain_B |...
            facelabel(:,1) == grain_B & facelabel(:,2) == grain_A);
        facetri_normal = tri_normal(mask_objface, :);
        angle_diffs = ones(size(facetri_normal, 1), 1)*90;

        % ----- Convert triangle normals to crystal frame and apply symmetries -----
%         g1 = reshape(G(face_to_calc(i, 1), :), 3,3)';
%         g2 = reshape(G(face_to_calc(i, 2), :), 3,3)';
        g1 = EAtoG(EA(grain_A, :));
        g2 = EAtoG(EA(grain_B, :));
        gg1 = O * g1;
        gg2 = O * g2;
        
        for j = 1:size(facetri_normal, 1)
            normal = facetri_normal(j, :)';
            % ##### Min Angle ##### 
            normal_g1 = reshape(gg1 * normal, 3, n_sym)';
            % ----- dist to (111) -----
            ang_diff_tmp1 = acosd(dot(normal_g1, fixed_normal_1, 2));
            if min(ang_diff_tmp1) < angle_diffs(j)
                angle_diffs(j) = min(ang_diff_tmp1);
            end 
            normal_g2 = reshape(gg2 * normal, 3, n_sym)';
            % ----- dist to (111) -----
            ang_diff_tmp2 = acosd(dot(normal_g2, fixed_normal_2, 2));
            if min(ang_diff_tmp2) < angle_diffs(j)
                angle_diffs(j) = min(ang_diff_tmp2);
            end
        end
        angle_diff_list = [angle_diff_list; angle_diffs];

    end
end




%% #######################################################################################  
% Test fixed normal 
g_fixed = AAToG(60, [1,1,1]);
fixed_normal_1 = [1,1,1]';
fixed_normal_1 = fixed_normal_1/norm(fixed_normal_1);
fixed_normal_2 = zeros(24*24, 3);

O = CrysSym;
idx = 1;
for i = 1:24
    for j = 1:24
        fixed_normal_2(idx, :) = O(:,:,i) * g_fixed * O(:,:,j) * fixed_normal_1;
        idx = idx + 1;
    end
end
        





