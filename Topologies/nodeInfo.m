function [numGBs, numTLs] = nodeInfo(objType, file, TLs, QNs, FCNs, SixCNs)
% ############################################################################
% * This function helps topology check by finding the #GBs and #TLs at a quadNode or superQuadNode.
% * Related files
%      - findTripleLines, findQuandNodes, findQuandNodes_checks
%      - Usually used in findQuandNodes_checks, better clean multiplicity QNs first.
% * objType = 'QNs' | 'FCNs' | 'SixCNs'
%      - Specifies the kind of nodes that needs to be find
% 
% * NOTE: QNs have to have been sorted rows
% ############################################################################

FLs = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
FLs(any(FLs<=0, 2), :) = [];
FLs = unique(FLs, 'rows');
FLs = sort(FLs, 2);

if strcmp(objType, 'QNs')
    % ----- #GBs at QN -----
    numGBs = zeros(length(QNs),1);
    C42 = combnk(1:4,2);
    for i = 1:length(QNs)
        sixGBs_QN = reshape(QNs(i, C42),6,2);
        numGBs(i) = sum(ismember(sixGBs_QN, FLs, 'rows'));
    end
    % ----- #bonds: #TLs at QN -----
    numTLs = zeros(length(QNs),1);
    C43 = combnk(1:4,3);
    for i = 1:length(QNs)
        fourTLs_QN = reshape(QNs(i, C43),4,3);
        numTLs(i) = sum(ismember(fourTLs_QN,TLs, 'rows'));
    end
    
elseif strcmp(objType, 'FCNs')
    % ----- #GBs at FNC -----
    numGBs = zeros(length(FCNs),1);
    C52 = combnk(1:5,2);
    for i = 1:length(FCNs)
        tenGBs_QN = reshape(FCNs(i, C52),10,2);
        numGBs(i) = sum(ismember(tenGBs_QN, FLs, 'rows'));
    end
    % ----- #bonds: #TLs at FCN -----
    numTLs = zeros(length(FCNs),1);
    C53 = combnk(1:5,3);
    for i = 1:length(FCNs)
        tenTLs_QN = reshape(FCNs(i, C53),10,3);
        numTLs(i) = sum(ismember(tenTLs_QN,TLs, 'rows'));
    end
    
elseif strcmp(objType, 'SixCNs')
    % ----- #GBs at FNC -----
    numGBs = zeros(length(SixCNs),1);
    C62 = combnk(1:6,2);
    for i = 1:length(SixCNs)
        tenGBs_QN = reshape(SixCNs(i, C62),15,2);
        numGBs(i) = sum(ismember(tenGBs_QN, FLs, 'rows'));
    end
    % ----- #bonds: #TLs at FCN -----
    numTLs = zeros(length(SixCNs),1);
    C63 = combnk(1:6,3);
    for i = 1:length(SixCNs)
        tenTLs_QN = reshape(SixCNs(i, C63),20,3);
        numTLs(i) = sum(ismember(tenTLs_QN,TLs, 'rows'));
    end
    
end

end
