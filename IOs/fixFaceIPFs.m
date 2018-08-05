function fixFaceIPFs(file)
% % ##################################################################
% Instruction
%   * This script is to fix the faceIPFcolor problem in current version of D3D. 
%       - The problem was due to D3D didn't have the winding of faceIPF
%       colors right. 
%       - To fix it, one simply revert the order of the
%       recorded IPFcolors. The correct order can be found to faceLabels.
%   * In order to run this script, one must have calculated the
%   faceIPFcolors from D3D and have named it IPFcolors. 
% ##################################################################
% ------------------------------------------------------
% Input for Debug
% file = '/Users/xiaotingzhong/Desktop/Datas/HCP twin/labeled_mesh2_forParaview.dream3d';
% ------------------------------------------------------
% ##### read in data #####
IPFs = h5read(file,'/DataContainers/TriangleDataContainer/FaceData/IPFColors').';
FLs = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';

% ##### find the indices for switching(idx_A) and not switching(idx_B) #####
FLs(any(FLs <= 0, 2), :) = 0;
mask_A = (IPFs(:,4)~=0 & IPFs(:,5)~=0 & IPFs(:,6)~=0 & FLs(:,1)>FLs(:,2));
mask_B = (IPFs(:,4)~=0 & IPFs(:,5)~=0 & IPFs(:,6)~=0 & FLs(:,1)<FLs(:,2));
idx = (1:size(IPFs, 1))';
idx_A = idx(mask_A);
idx_B = idx(mask_B);

% ##### switch IPF color order #####
holder = zeros(size(IPFs, 1), 3);
holder(idx_A, :) = IPFs(idx_A, 1:3);
holder(idx_B, :) = IPFs(idx_B, 4:6);

IPFs(idx_A, 1:3) = IPFs(idx_A, 4:6);
IPFs(:, 4:6) = holder;

% ##### write back to D3D file #####
h5write(file,'/DataContainers/TriangleDataContainer/FaceData/IPFColors', IPFs');

end