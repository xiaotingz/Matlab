tmp = h5read('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth_croptrial.dream3d', ...
    '/DataContainers/ImageDataContainer/CellData/FeatureIds');
data = squeeze(tmp);

for i = 1:size(data, 3)
    fname = ['An5new6_smooth_grainID_', num2str(i), '.raw'];
    fid=fopen(fname,'w+');
    cnt=fwrite(fid,data(:,:,i), 'uint16');
    fclose(fid);
end


