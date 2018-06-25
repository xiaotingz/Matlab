function [] = smallQNInfo(idx, smallQN, QNs, TLs)
% ############################################################################
% - This function displays ID of the QuadNode and IDs of the associated tripleLines
% - Used in findQuadNodes_checks
% ############################################################################
C43 = combnk(1:4,3);

idx_SQN = smallQN(idx);
SQN_IDs = QNs(idx_SQN, :);
SQN_TLs = reshape(QNs(idx_SQN, C43),4,3);
SQN_TLs = SQN_TLs(ismember(SQN_TLs, TLs, 'rows'), :);
fprintf(['\n quadNodeID: \n', num2str(SQN_IDs), '\n']);
fprintf('\n associated TLs: \n');
for i = 1:size(SQN_TLs, 1)
    disp(num2str(SQN_TLs(i, :)));
end
fprintf('\n');

end