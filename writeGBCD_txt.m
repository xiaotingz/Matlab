% SampleGBCD = textread('SampleGBCD(Ni_Bi_gbcd_shifted_skip).txt');
% GBCD_raw = h5read('/Volumes/RESEARCH/Simple Geometry Jun.21_Aug/Aug.24 r1_v4abs/Aug24_r1_GBCD.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
% 
% GBCD = GBCD_raw(:,1);
% 
% fileID = fopen('Aug24_r1_GBCD.txt','w');
% fprintf(fileID,'%12.8f\n',GBCD);
% fclose(fileID);

GBCD = textread('Aug24_r1_GBCD.txt');
 
GBCD_flip = flipud(GBCD);
fileID = fopen('Aug24_r1_GBCD_flip.txt','w');
fprintf(fileID,'%12.8f\n',GBCD_flip);
fclose(fileID);
% 
% GBCD_flip = textread('Aug24_r1_GBCD_flip.txt');
% 
% j = 1;
% for i = 1:length(GBCD)
%     if GBCD(i) ~= 0
%         gbcd(j) = GBCD(i);
%         j = j + 1;
%     end
% end
% j = 1;
% for i = 1:length(GBCD_flip)
%     if GBCD_flip(i) ~= 0
%         gbcd_flip(j) = GBCD_flip(i);
%         j = j + 1;
%     end
% end