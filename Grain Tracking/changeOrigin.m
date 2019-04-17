file_An4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6.dream3d');
file_An5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh.dream3d');

origin_An4 = double(h5read(file_An4,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN'));
origin_An5 = double(h5read(file_An5,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN'));

origin_An4_new = origin_An5;
h5write('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_changeOrigin.dream3d','/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', origin_An4_new);







