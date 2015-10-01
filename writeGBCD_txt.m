GBCD_raw = h5read('/Volumes/RESEARCH/Simple Geometry Jun.21_Aug/Aug.24 s4_v4abs/Aug24_s4_GBCD.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');

GBCD = GBCD_raw(:,1);

fileID = fopen('s4_GBCDRes9.txt','w');
fprintf(fileID,'%12.8f\n',GBCD);
fclose(fileID);

