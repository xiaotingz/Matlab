% load('/Volumes/XIAOTING/Ni/working/190621_tracked_faces_full.mat', 'face_area_an4', ...
%     'face_area_diff', 'face_itg_abscurv_an4', 'face_itg_abscurv_diff','tracked_uniqueface_an4_full');
% file = file_an4;
face_centroids = calcFaceCentroid(file_an4, tracked_uniqueface_an4);

% load('190627_an4an5_Lsmooth_tmp.mat')
face_area_diff_ratio = face_area_diff ./ face_area_an4;
mask = face_area_an4 > 20 & (face_area_an4 + face_area_diff) > 20 & ...
    face_area_diff_ratio < 10 & face_area_diff_ratio > -0.9;
face_area_an4 = face_area_an4(mask);
face_area_diff = face_area_diff(mask);
face_centroids = face_centroids(mask, :);
face_itg_abscurv_an4 = face_itg_abscurv_an4(mask);
face_itg_abscurv_diff = face_itg_abscurv_diff(mask);
tracked_uniqueface_an4 = tracked_uniqueface_an4(mask, :);
%% -------------- Area & ?Area --------------
figure('Renderer', 'painters', 'Position', [10 10 1000 400])
ax1 = subplot(1, 2, 1);
x = face_area_an4_full;
coords = face_centroids_full;
thres = prctile(x,95);
mask = x > thres;
x = x(mask);
coords = coords(mask, :);
scale = (x - min(x))/max(x)*500 + 0.1;
scatter3(coords(:, 1), coords(:, 2), coords(:, 3), ...
    scale, scale, 'filled');
daspect([1 1 1])
title('top 5% area', 'fontSize', 20)
cb1 = colorbar('westoutside');

ax2 = subplot(1, 2, 2);
x = face_area_diff_full;
coords = face_centroids_full;
thres = prctile(x,95);
mask = x > thres;
x = x(mask);
coords = coords(mask, :);
scale = (x - min(x))/max(x)*500 + 0.1;
scatter3(coords(:, 1), coords(:, 2), coords(:, 3), ...
    scale, scale, 'filled');
daspect([1 1 1])
title('top 5% area change', 'fontSize', 20)
cb2 = colorbar('eastoutside');

% colormap('jet')
rotate3d on
Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);


%% -------------- (Face Integral Curvature) & ?(Face Integral Curvature) --------------
figure('Renderer', 'painters', 'Position', [10 10 1000 400])
ax1 = subplot(1, 2, 1);
x = face_itg_abscurv_an4_full;
coords = face_centroids_full;
thres = prctile(x,95);
mask = x > thres;
x = x(mask);
coords = coords(mask, :);
scale = (x - min(x))/max(x)*500 + 0.1;
scatter3(coords(:, 1), coords(:, 2), coords(:, 3), ...
    scale, scale, 'filled');
daspect([1 1 1])
title('top 5% integral curvature', 'fontSize', 20)
cb1 = colorbar('westoutside');

ax2 = subplot(1, 2, 2);
x = face_itg_abscurv_diff_full;
coords = face_centroids_full;
thres = prctile(x,95);
mask = x > thres;
x = x(mask);
coords = coords(mask, :);
scale = (x - min(x))/max(x)*500 + 0.1;
scatter3(coords(:, 1), coords(:, 2), coords(:, 3), ...
    scale, scale, 'filled');
daspect([1 1 1])
title('top 5% integral curvature change', 'fontSize', 20)
cb2 = colorbar('eastoutside');

% colormap('jet')
% colorbar
rotate3d on
Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);






