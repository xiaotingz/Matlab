for j = 38:78
    DigitTen = floor(j/10);
    DigitOne = rem(j,10);
    InFile = strcat('/Users/xiaotingzhong/Desktop/STO_1470/subset2/SrTiO3_0', int2str(DigitTen),int2str(DigitOne),'_Mod.ang');
    data = textread(InFile);
    
    OutFile = strcat('/Users/xiaotingzhong/Desktop/STO_1470/subset1_Modified/SrTiO3_0', int2str(DigitTen),int2str(DigitOne),'.ang');
    fileID = fopen(OutFile,'w');
    fprintf(fileID,'# TEM_PIXperUM          1.000000\n');
    fprintf(fileID,'# x-star                0.450410\n');
    fprintf(fileID,'# y-star                0.799933\n');
    fprintf(fileID,'# z-star                0.699186\n');
    fprintf(fileID,'# WorkingDistance       5.000000\n');
    fprintf(fileID,'#\n');
    fprintf(fileID,'# Phase 1\n');
    fprintf(fileID,'# MaterialName  	Strontium Titanium Oxide\n');
    fprintf(fileID,'# Formula     	SrTiO3\n');
    fprintf(fileID,'# Info \n');
    fprintf(fileID,'# Symmetry              43\n');
    fprintf(fileID,'# LatticeConstants      3.905 3.905 3.905  90.000  90.000  90.000\n');
    fprintf(fileID,'# NumberFamilies        5\n');
    fprintf(fileID,'# hklFamilies   	 1  1  0 1 0.000000 1\n');
    fprintf(fileID,'# hklFamilies   	 1  1  1 1 0.000000 1\n');
    fprintf(fileID,'# hklFamilies   	 2  0  0 1 0.000000 1\n');
    fprintf(fileID,'# hklFamilies   	 2  1  1 1 0.000000 1\n');
    fprintf(fileID,'# hklFamilies   	 3  1  0 1 0.000000 1\n');
    fprintf(fileID,'# Categories16993468 -1275068416 2009465508 2009385371 0\n');
    fprintf(fileID,'#\n');
    fprintf(fileID,'# GRID: SqrGrid\n');
    fprintf(fileID,'# XSTEP: 0.300000\n');
    fprintf(fileID,'# YSTEP: 0.300000\n');
%     CHANGE for different subsets
    fprintf(fileID,'# NCOLS_ODD: 225\n');
    fprintf(fileID,'# NCOLS_EVEN: 225\n');
    fprintf(fileID,'# NROWS: 312\n');
    fprintf(fileID,'#\n');
    fprintf(fileID,'# OPERATOR: 	supervisor\n');
    fprintf(fileID,'#\n');
    fprintf(fileID,'# SAMPLEID: \n');
    fprintf(fileID,'#\n');
    fprintf(fileID,'# SCANID:\n');
    fprintf(fileID,'#\n');

    format = '%8.3f    %8.3f    %8.3f    %8.3f    %8.3f  %6.2f  %6.3f  %d  %d\n';
    for i = 1:length(data)
        if data(i,7) < 0
            data(i,7) = 0;
        end
        fprintf(fileID,format,data(i,1), data(i,2), data(i,3), data(i,4),data(i,5),data(i,6),data(i,7),data(i,8),data(i,9));
    end

end


