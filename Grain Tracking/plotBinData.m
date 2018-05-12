function plotBinData(x, y, xrange, yrange, stepsize, label_x, label_y, backgroundOption, STDoption)
% ###########################################################################
% * misA UNIT: DEGREE!!
% * Output
%     - data_grid = [leftBoundaryOfCurrentBin, avg(y), countInBin]
% #########################################################################
% ------------------ load data for debug --------------------
% clear
% x = plotData_6(:,1);
% y  = plotData_6(:,2);
% xrange = [-30, 20];
% % yrange = [-1e4, 1e4];
% stepsize = 1;
% label_x = 'Face #Edges Difference';
% label_y = 'Face Area Difference';
% STDoption = false;
% backgroundOption = true;
% showCnt = backgroundOption;
% showSTD = STDoption;
% x = 1:1 : 50;
% xrange = [1, 50];
% y = 1:1:50;
% stepsize = 3;
% -----------------------------------------------------------
x = double(x);
y = double(y);

if nargin > 8
    showCnt = backgroundOption;
    showSTD = STDoption;
elseif nargin > 7 && nargin < 8
    showCnt = backgroundOption;
    showSTD = true;
else
    showCnt = false;
    showSTD = true;
end

mask = ((x > xrange(1) & x < xrange(2)) & (y > yrange(1) & y < yrange(2)));
x = x(mask);
y = y(mask);

% if xrange(1) < 0
%     xPlot = [-flip(0:stepsize:xrange(2)),stepsize:stepsize:xrange(2)];
% else
xPlot = [xrange(1):stepsize:xrange(2)];
% end

data_grid = zeros(length(xPlot), 5);
data_grid(:,1) = xPlot;
for i = 1:length(data_grid)-1
    lowLim = data_grid(i);
    upLim = data_grid(i + 1);
    mask2 = (x >= lowLim & x < upLim);
    data_grid(i,2) = sum(x(mask2))/sum(mask2);
    data_grid(i,3) = sum(y(mask2))/sum(mask2);
    data_grid(i,4) = sum(mask2);
    if isempty(y(mask2))
        data_grid(i,5) = 0;
    else
        data_grid(i,5) = std(y(mask2));
    end
    if showSTD 
        line([data_grid(i,2), data_grid(i,2)], [data_grid(i,3)-data_grid(i,5), data_grid(i,3)+data_grid(i,5)], 'Color','k')
        hold on
    end
end
% -- if there is no Y data, there is no data in the bin -- 
mask_nodata = (data_grid(:,3) == 0);
data_grid(mask_nodata,:) = [];
scatter(data_grid(:,2), data_grid(:,3), 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
line([xrange(1), xrange(2)], [0,0], 'LineStyle', '--', 'Color', [0.5, 0.5, 0.5]);
% line([0,0], [-1500, 3000], 'LineStyle', '--', 'Color', [0.5, 0.5, 0.5]);
xlabel(label_x);
ylabel(label_y);
% ylim([-4000,6000]);
box on

if showCnt
    yyaxis right
%     histogram('BinEdges', [data_grid(1,1)-stepsize; data_grid(:,1)], 'BinCounts',data_grid(:,4))
    bar(data_grid(:,2), data_grid(:,4),'Barwidth', 1.1,'FaceColor', [0.5,0.5,0.5], 'EdgeColor', [0.5,0.5,0.5], 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
%     ylabel('Counts')
    ylabel('#Faces in Bin')
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    xlim(xrange);
end

end
