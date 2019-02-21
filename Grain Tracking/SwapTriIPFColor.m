file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth_forParaview.dream3d';
IPF_switch = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/IPFColors')).';
facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';


for i = 1:length(IPF_switch)
    if facelabel(i,1) < facelabel(i,2)
        tmp = IPF_switch(i, 4:6);
        IPF_switch(i, 4:6) = IPF_switch(i, 1:3);
        IPF_switch(i, 1:3) = tmp;
    end
end


h5write(file,'/DataContainers/TriangleDataContainer/FaceData/IPFColors_BlueFirst', IPF_switch');
h5write(file,'/DataContainers/TriangleDataContainer/FaceData/IPFcolors_Blue', IPF_BlueDiff');


%%
file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
istwin = boolean(h5read(file, '/DataContainers/TriangleDataContainer/FaceData/TwinBoundary'));
min_angles_all = 180*ones(size(istwin));

min_angles_all(istwin) = min_angles;

h5write('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth_forParaview.dream3d'...
    ,'/DataContainers/TriangleDataContainer/FaceData/TwinBoundary', int8(istwin));
h5write('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth_forParaview.dream3d'...
    ,'/DataContainers/TriangleDataContainer/FaceData/CTwin', min_angles_all);







