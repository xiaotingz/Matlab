%% ##### PLOT ##### 
% ----- EC = [(1:length(F))', F, E, V, V-E+F] ----- 

histogram(EC_An4(:,4))
hold on
histogram(EC_An5(:,5))


xlabel(label_x,'Interpreter','latex');
ylabel(label_y,'Interpreter','latex');



%% ##### check CI and change voxel orientation #####
file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4_cropSixCN_886184.dream3d';
[result, QNList, fiveCoordNList] = findQuadNodes(file)

