set(0,'defaultAxesFontSize',18)

% --------------------------- load data ---------------------------
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/final/190730_Hsmoth_GeoTopo.mat',...
    'face_area_an4', 'face_area_diff', 'face_itg_abscurv_diff')
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/final/190730_Hsmoth_Other.mat', ...
    'dists_miso_an4', 'dists_norm_an4', 'std_norm_an4')
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/final/190730_Hsmoth_MeanField.mat', ...
    'tracked_uniqueface_an4_complete', 'tracked_uniqueface_an4_inner')
% --------------------------- clean data ---------------------------
face_area_an5 = face_area_an4 + face_area_diff;
face_area_diff_ratio = face_area_diff ./ face_area_an4;
mask_complete = ismember(tracked_uniqueface_an4_inner, tracked_uniqueface_an4_complete, 'rows');
mask_good = face_area_an4 > 20 & face_area_an5 > 20 & face_area_diff_ratio < 10 & face_area_diff_ratio > -0.9;
mask = mask_complete & mask_good;

face_area_an4 = face_area_an4(mask);
face_area_diff = face_area_diff(mask);
face_itg_abscurv_diff = face_itg_abscurv_diff(mask);
dists_miso_an4 = dists_miso_an4(mask, :);
% dists_miso_an5 = dists_miso_an5(mask, :);
dists = dists_miso_an4;

% --------------------------- box plot ---------------------------
data = [face_area_diff, face_itg_abscurv_diff];
data_sigma3 = data(dists(:,1) < 5, :);
data_sigma5 = data(dists(:,2) < 5, :);
data_sigma7 = data(dists(:,3) < 5, :);
data_sigma9 = data(dists(:,4) < 5, :);
data_sigma27a = data(dists(:,5) < 5, :);
data_sigma27b = data(dists(:,6) < 5, :);
data_general = data(all(dists>5, 2), :);

% labels = {['#(sigma3) = ', num2str(length(data_sigma3))], ['#(sigma5) = ', num2str(length(data_sigma5))], ...
%         ['#(sigma7) = ', num2str(length(data_sigma7))], ['#(sigma9) = ', num2str(length(data_sigma9))], ...
%         ['#(sigma27a) = ', num2str(length(data_sigma27a))], ['#(sigma27b) = ', num2str(length(data_sigma27b))],...
%         ['#general = ', num2str(length(data_general))]};

miso_names = {'$\Sigma3$','$\Sigma5$','$\Sigma7$','$\Sigma9$','$\Sigma27a$','$\Sigma27b$', 'Others'};
miso_populations = sum(dists < 5);
miso_populations = [miso_populations, length(data) - sum(miso_populations)];

idx_plot = 1;
area_diff = [data_sigma3(:,idx_plot); data_sigma5(:,idx_plot); data_sigma7(:,idx_plot); data_sigma9(:,idx_plot); ...
            data_sigma27a(:,idx_plot); data_sigma27b(:,idx_plot); data_general(:,idx_plot)];
miso_id = [ones(size(data_sigma3(:,1))); ones(size(data_sigma5(:,1)))*2; ...
            ones(size(data_sigma7(:,1)))*3; ones(size(data_sigma9(:,1)))*4; ...
            ones(size(data_sigma27a(:,1)))*5; ones(size(data_sigma27b(:,1)))*6; ...
            ones(size(data_general(:,1)))*7];
        
% figure('Renderer', 'painters', 'Position', [10 10 900 500])
% ----- Plot -----
boxplot(area_diff, miso_id, 'Symbol', '', 'Colors', 'k');
hold on
line ([0, 8], [0, 0], 'Color', [0.5, 0.5, 0.5], 'LineStyle',':')
ylim([-650, 650])

% ----- Set double x-axis -----
ax1 = gca;
set(ax1,'XTickLabel', miso_names,'TickLabelInterpreter', 'latex');
% xtickangle(45)
ylabel('$\Delta$A ($\mu m ^{2}$)', 'Interpreter', 'latex')
xlabel('Misorientation', 'Interpreter', 'latex')
ax2 = axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation','right','Color','none');
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
set(ax2,'XTickLabel', miso_populations,'YTickLabel', [],'TickLabelInterpreter', 'latex');
xlabel('Population', 'Interpreter', 'latex')
box on
% ----- Adjust plot field of view -----
TightInset = max([ax1.TightInset; ax2.TightInset]); %[Left, Bot, Right, Top] border needed
AxPos = [TightInset(1)+0.01, TightInset(2), 1-sum(TightInset([1 3]))-0.02, 1-sum(TightInset([2 4]))-0.02]; %New position
ax1.Position = AxPos;
ax2.Position = AxPos;
% ----- Save figure -----
print('miso_distribution','-dtiff','-r300')



% disp([num2str(sum(face_area_an4)), ',  ', num2str(sum(face_area_an5)), ',  ', ...
%     num2str((sum(face_area_an5) - sum(face_area_an4)) / sum(face_area_an4))])
% disp([num2str((sum(face_area_an4(dists(:,1) < 5)))), ',  ', num2str((sum(face_area_an5(dists(:,1) < 5)))), ',  ', ...
%     num2str((sum(face_area_an5(dists(:,1) < 5)) - sum(face_area_an4(dists(:,1) < 5))) / sum(face_area_an4(dists(:,1) < 5)))])
%%
% --------------------------- cdf plot ---------------------------
figure()
h(1) = cdfplot(data_sigma3(:,1));
hold on
h(2) = cdfplot(data_sigma5(:,1));
h(3) = cdfplot(data_sigma7(:,1));
h(4) = cdfplot(data_sigma9(:,1));
h(5) = cdfplot(data_sigma27a(:,1));
h(6) = cdfplot(data_sigma27b(:,1));
h(7) = cdfplot(data_general(:,1));
set( h(:), 'LineWidth',2);
xlim([-1000, 1000])
legend(labels, 'Location', 'southeast')
xlabel('dA')
% title('all tracked faces')
% title('no small faces')
title('full tracked faces')
print('dA_byMisorientation_cdf_full_tracked','-dtiff','-r300')

% csvwrite('no_small_face_mask.txt', mask)

