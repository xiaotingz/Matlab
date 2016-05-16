GBCD_raw = h5read('/Users/xiaotingzhong/Desktop/Datas/setTo0/Jan31_Aca0_CurvDistri_10.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');

GBCD = GBCD_raw(:,1);

fileID = fopen('A_CurvDistri10.txt','w');
fprintf(fileID,'%12.8f\n',GBCD);
fclose(fileID);

