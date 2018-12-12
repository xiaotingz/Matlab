file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth_forParaview.dream3d';
IPF_switch = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/IPFColors_BlueFirst')).';
IPF_maxBlue = zeros(size(IPF_switch,1), 1);

for i = 1:length(IPF_switch)
    if IPF_switch(i, 6) > IPF_switch(i, 3)
        tmp = IPF_switch(i, 1:3);
        IPF_switch(i, 4:6) = IPF_switch(i, 1:3);
        IPF_switch(i, 1:3) = tmp;
    end
    IPF_maxBlue(i) = IPF_switch(i, 3);
end


%%
h5write(file,'/DataContainers/TriangleDataContainer/FaceData/IPFColors_BlueFirst', int8(IPF_switch'));
h5write(file,'/DataContainers/TriangleDataContainer/FaceData/IPFcolors_Blue', int8(IPF_maxBlue'));