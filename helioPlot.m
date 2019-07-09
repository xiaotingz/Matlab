% ##### normalize the coefficients #####
% """
% The base circle is 0.6, with bars of heights +- 0.4
% """
% ----- read data -----
file_name = fopen('/Users/xiaotingzhong/Documents/Matlab/Grain Tracking/_DfMs_XGB_feature_names.txt','r');
file_coef = fopen('/Users/xiaotingzhong/Documents/Matlab/Grain Tracking/_DfMs_XGB_coefs.txt','r');
fig_name = 'DfMs_XGB';
% feature_class_idx = [[1, 6]; [7, 8]; [9, length(feature_names)-6]; [length(feature_names)-5, length(feature_names)]];
% feature_class_idx = [[1, 4]; [5, 6]; [7, length(feature_names)-6]; [length(feature_names)-5, length(feature_names)]];
feature_class_idx = [[1, 4]; [5, 8]; [9, 18]; [19, 21]];

coefs = textscan(file_coef,'%f');
coefs = coefs{1,1};
feature_names = textscan(file_name,'%s');
feature_names = feature_names{1,1};

n = length(coefs);
base = 0.6;
coef_scaled = coefs'./max(abs(coefs)) * (1-base);
head = base*ones(n, 1)';
tail = head + coef_scaled;

% ----- adjust figure range -----
% """ plot an invisible figure to make everything in plot """
h = figure(1);
set(h, 'Position', [100, 100, 800, 800])
polarplot([0; 0], [0; 1.8], 'color', 'w')
hold on

% ----- plot the bars -----
theta = deg2rad(0:360/(n):360-0.1);
polarplot([theta; theta], [head; tail], 'lineWidth',12, 'color', 'k')

% ----- change the color of negative bars -----
edge = 0.003;
mask = coef_scaled < 0;
polarplot([theta(mask); theta(mask)], [head(mask)-edge; tail(mask)+edge], 'lineWidth',11, 'color', 'w')

% ----- base circle -----
grid off
pax = gca;
pax.ThetaTickLabel = [];
pax.RTickLabel = [];
base_angs = linspace(0,2*pi,360);
base_r = base*ones(size(base_angs));
polarplot(base_angs, base_r, 'color', [0.5, 0.5, 0.5])

% ----- labels -----
r_out = 1.15;
for i = 1:length(theta)
    ang = rad2deg(theta(i));
    if ang > 90 && ang < 270
        ang = ang + 180;
        text(theta(i), r_out,  char(feature_names{i, 1}),'Rotation', ang, ...
            'HorizontalAlignment', 'right', 'fontSize', 12, 'Interpreter', 'none')
    else
        text(theta(i), r_out,  char(feature_names{i, 1}), 'Rotation', ang, ...
            'fontSize', 12, 'Interpreter', 'none')
    end
end

% ----- labels lines -----
% out_r = out_r*ones(size(base_angs));
polarplot([theta; theta], [base.*ones(size(theta)); r_out.*ones(size(theta))], ':', 'color', [0.5, 0.5, 0.5])

% ----- feature class: ticks and archs -----
r_out_2 = r_out - 0.05;
r_out_3 = r_out - 0.02;
for i = 1:size(feature_class_idx)
    ang_1 = theta(feature_class_idx(i, 1));
    ang_2 = theta(feature_class_idx(i, 2));
    polarplot([ang_1; ang_1], [r_out_2; r_out_3], 'lineWidth', 1,'color', 'k');
    polarplot([ang_2; ang_2], [r_out_2; r_out_3], 'lineWidth', 1,'color', 'k');
    base_angs = linspace(ang_1, ang_2, round(rad2deg(ang_2 - ang_1)));
    polarplot(base_angs, r_out_2*ones(size(base_angs)), 'lineWidth', 1, 'color', 'k');
end

% ----- save_plot -----
print(fig_name,'-dtiff','-r300')



