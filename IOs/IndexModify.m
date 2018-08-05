% First, get the index for the two system.
    % phi_1,PHI,phi_2,theta,Phi
    % D3D is increasing from right to left, gbcd_graph is increasing from left to right.
NumStep = 10;
index = zeros(NumStep^5,5);

    % list all the bins as the order in gbcd_graph. 
    % at the same time, calculate gbcd_index using rule in D3D, namely order in D3D.
n1 = NumStep;
n1n2 = NumStep^2;
n1n2n3 = NumStep^3;
n1n2n3n4 = NumStep^4;
i = 1;
    % NOTICE, Phi is from 0 to 2Pi,so the size of e is 4*NumStep 
for a = 1 : NumStep
    for b = 1 : NumStep
        for c = 1 : NumStep
            for d = 1 : NumStep
                for e  = 1 : NumStep*4
                    index(i,1) = a;
                    index(i,2) = b;
                    index(i,3) = c;
                    index(i,4) = d;
                    index(i,5) = e;
                    
%                     calculate the index in D3D system -- gbcd_index 
                    gbcd_index(i) = a + b*n1 + c*n1n2 + d*n1n2n3 + e*n1n2n3n4;
                    i = i + 1;
                end
            end
        end
    end
end

    % matlab begins with 1 so there are side effects
% 7380 for step=10 degree; 11110 for step=9 degree
gbcd_index = (gbcd_index - 11110).';
    % the index matrix is in the order of graph_gbcd, so the corresponding graph_index is the natural order.
graph_index = unique(gbcd_index);
%%
% load the index file according to the GBCD resolution used
load('TwoIndex_res10.mat')
    % link the two system. Sort as gbcd_index in its natural order, then graph_index can be used to place D3D gbcd data
d3d_to_graph = [gbcd_index,graph_index];
d3d_to_graph = sortrows(d3d_to_graph);


% % ######################## Write GBCD ########################
% data_read = textread('/Users/xiaotingzhong/Desktop/programs/energy/725_project/iter300/energy_v.res');
% GBCD_raw = h5read('/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/Jan31_Fca0_gbcd10.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
% data_read = GBCD_raw(:,1);
% 
% data_converted = zeros(length(data_read),1);
% for i = 1:length(data_read)
%     data_converted(d3d_to_graph(i,2)) = data_read(i);
% end
% 
% fileID = fopen('/Users/xiaotingzhong/Desktop/F_1volume_GBCD.txt','w');
% fprintf(fileID,'%12.8f\n',data_converted);
% fclose(fileID);


% % ######################## Write GBCurvD ########################
GBCurvD_1 = h5read('/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/problemData/Jan31_Aa0_CurvDistri10.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
GBCurvD_2 = h5read('/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/problemData/Jan31_Aa0_CurvDistri10.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCDCounters');
data_read = [GBCurvD_1(:,1), GBCurvD_2(:,1)];
data_converted = zeros(length(data_read),2);
for i = 1:length(data_read)
    data_converted(d3d_to_graph(i,2), :) = data_read(i, :);
end

fileID = fopen('/Users/xiaotingzhong/Desktop/A_GBCurvD_Jan31_Aa0_CurvDistri10.txt','w');
for i = 1:length(data_converted)
    fprintf(fileID,'%12.8f %12.8f\n',data_converted(i,1), data_converted(i,2));
end
fclose(fileID);

