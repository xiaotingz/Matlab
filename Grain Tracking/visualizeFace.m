% file_coordsmooth_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
% file_coordsmooth_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
% file_1 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_forParaview.dream3d');
% file_2 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh_forParaview.dream3d');
obj_face = face_node_info_3;

figure(1)
trisurf(reshape(obj_face{2,1}, 3, [])', node_coordmesh_an4(:,1), node_coordmesh_an4(:,2), node_coordmesh_an4(:,3),'Facecolor',[0, 0.4470, 0.7410] );
hold on
trisurf(reshape(obj_face{2,2}, 3, [])', node_coordmesh_an5(:,1), node_coordmesh_an5(:,2), node_coordmesh_an5(:,3),'Facecolor',[0.9290, 0.6940, 0.1250]);
rotate3d on

figure(2)
trisurf(reshape(obj_face{2,1}, 3, [])', node_coordsmooth_an4(:,1), node_coordsmooth_an4(:,2), node_coordsmooth_an4(:,3),'Facecolor',[0, 0.4470, 0.7410] );
hold on
trisurf(reshape(obj_face{2,2}, 3, [])', node_coordsmooth_an5(:,1), node_coordsmooth_an5(:,2), node_coordsmooth_an5(:,3),'Facecolor',[0.9290, 0.6940, 0.1250]);
rotate3d on