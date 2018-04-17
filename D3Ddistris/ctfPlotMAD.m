file = ('/Users/xiaotingzhong/Desktop/Datas/Ti/step05_mesh.dream3d');
MAD = double(h5read(file,'/DataContainers/ImageDataContainer/CellData/EulerAngles'));


h = histogram(MAD(1,:,:,:));
h.Normalization = 'probability';
h.NumBins=100;
xlim([0,3.5])

ax = gca;
set(ax,'fontsize',19)
xlabel('MAD','FontSize',21);
ylabel('probability','FontSize',21);
