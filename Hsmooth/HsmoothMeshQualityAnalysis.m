
% set(0,'defaultAxesFontSize',20)

% file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin2_smooth.dream3d';
% file_an5 = '/Volumes/XIAOTING/Ni/An5new6_smooth.dream3d';
file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
file_an5 = '/Volumes/XIAOTING/Ni/An5new6_Hsmooth.dream3d';
fl_an4 = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels');
fl_an5 = h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels');
min_ang_an4 = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
min_ang_an5 = roundn(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
area_an4 = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
area_an5 = roundn(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
curv_an4 = abs(roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5))';
curv_an5 = abs(roundn(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5))';

% mask_an4 = all(fl_an4>0, 2);
% min_ang_an4 = min_ang_an4(mask_an4);
% area_an4 = area_an4(mask_an4);
% curv_an4 = curv_an4(mask_an4);
% mask_an5 = all(fl_an5>0, 2);
% min_ang_an5 = min_ang_an5(mask_an5);
% area_an5 = area_an5(mask_an5);
% curv_an5 = curv_an5(mask_an5);

eps_ang = 10;
eps_curv = 1;
eps_area = 7;
mask_an4 = (min_ang_an4 > eps_ang & curv_an4 < eps_curv);
mask_an4 = (mask_an4 & area_an4 < eps_area);
mask_an5 = (min_ang_an5 > eps_ang & curv_an5 < eps_curv);
mask_an5 = (mask_an5 & area_an5 < eps_area);
good_tris_an4 = double(sum(mask_an4)) / double(length(min_ang_an4));
good_tris_an5 = double(sum(mask_an5)) / double(length(min_ang_an5));
disp(['good_triangles in an4 = ', num2str(good_tris_an4), ',  good_triangles in an5 = ', num2str(good_tris_an5)])



%%
plot_an4 = min_ang_an4;
plot_an5 = min_ang_an5;
fname = 'Hsmooth triMinAngDist';
x_axis = 'min(mesh triangle inner angle)';


% histogram(area_an4(area_an4<10))
% hold on
% histogram(area_an5(area_an5<10))
% 
% legend('an4', 'an5')

h1 = cdfplot(area_an4);
hold on
h2 = cdfplot(area_an5);

set( h1,'LineWidth',3);
set( h2,'LineWidth',3, 'LineStyle', '--');
legend(['an4, 1 percentile = ', num2str(prctile(plot_an4, 1))], ...
    ['an5, 1 percentile = ', num2str(prctile(plot_an5, 1))], ...
    'Location','northwest'); 
%     'Location','southeast');
    
xlabel(x_axis)



title(fname)
print(fname, '-dtiff', '-r300')