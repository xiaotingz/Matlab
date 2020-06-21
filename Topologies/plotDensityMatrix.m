function plotDensityMatrix(x, y, edges_x, edges_y, x_name, y_name)
% ############################################################################
% * Inputs
%   - x & y
%       The two data arrays. Each row should be in correspondence. 
% * Note
%   - This function takes in two data arrays and plot their density matrix.
%   - See /Grain Curvature/D_F_plot_Kris.m for original code.
% ############################################################################
% ----------------------- load debug data -----------------------
% x = num_corners_an4;
% y = num_edges_an4;
% xbins = 35; xMin = 0; xMax = 35;
% ybins = 35; yMin = 0; yMax = 35;
% ---------------------------------------------------------------
% xMin = x_bin_setting(1);
% xMax = x_bin_setting(2);
% xbins = x_bin_setting(3);
% yMin = y_bin_setting(1);
% yMax = y_bin_setting(2);
% ybins = y_bin_setting(3);

% ----- Make Plot -----
% N = hist3([y,x], 'edges', {yMin:(yMax-yMin)/ybins:yMax xMin:(xMax-xMin)/xbins:xMax} );
N = hist3([y,x], 'edges', {edges_y edges_x} );
colormap(flipud(hot))
set(gca,'YDir','normal');
p = imagesc(N);
hold on

% ----- Decorate Edges -----
for i = 1:size(N,1)
    for j = 1:size(N,2)
        if N(i,j) ~= 0
            line([j-0.5, j-0.5],[i-0.5, i+0.5],'color', [0.8, 0.8, 0.8])
            line([j+0.5, j+0.5],[i-0.5, i+0.5],'color', [0.8, 0.8, 0.8])
            line([j-0.5, j+0.5],[i-0.5, i-0.5],'color', [0.8, 0.8, 0.8])
            line([j-0.5, j+0.5],[i+0.5, i+0.5],'color', [0.8, 0.8, 0.8])
            hold on
        end
    end
end
% yTicks = [0 : 10*ybins/55 : ybins*50/55];
% xTicks = [0, xbins/9: 2*xbins/9: 50];
% yLabels = ({'0', '10', '20', '30', '40', '50'});
% xLabels = ({'1', '2', '4', '6', '8', '10'});
% set(gca, 'YTick', yTicks, 'XTick', xTicks, 'XTickLabels',xLabels, 'YTickLabels',yLabels,  'XMinorTick','on','YMinorTick','on')

% ----- Finish settings -----
pbaspect([1 1 1])
axis('xy')
xlabel(x_name,'FontSize',20);
ylabel(y_name,'FontSize',20);
set(gca,'fontsize',19)
colorbar

end