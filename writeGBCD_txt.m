GBCD_raw = h5read('/Volumes/RESEARCH/Simple Geometry Jun.21_Aug/Aug.24 r1_v4abs/Aug24_r1_CurvDistri.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');

GBCD = GBCD_raw(:,1);

fileID = fopen('r1_CurvDistriRes10.txt','w');
fprintf(fileID,'%12.8f\n',GBCD);
fclose(fileID);

