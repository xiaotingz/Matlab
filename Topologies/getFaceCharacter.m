function [numCorners, numEdges] = getFaceCharacter(labels, TLs, QNs, FCNs)
% ############################################################################
% - labels = [n, 2]
% - This function takes only FiveCoordSuperQN but not other superQNs.
% ############################################################################
% ----------------------- load debug data -----------------------
% A = [277, 277, 277]';
% B = [8, 69, 100]';
% labels = [277, 8];
% ---------------------------------------------------------------
numEdges = zeros(size(labels, 1), 1);
numCorners = zeros(size(labels, 1), 1);

% ##### if only QuadPoints #####
if nargin == 3
    for i = 1:size(labels, 1)
        A = labels(i, 1);
        B = labels(i, 2);

        mask_TL_A = (TLs == A);
        mask_TL_B = (TLs == B);
        mask_TL_AB = ((sum(mask_TL_A, 2) + sum(mask_TL_B, 2)) == 2);
        numEdges(i) = sum(mask_TL_AB);

        mask_QN_A = (QNs == A);
        mask_QN_B = (QNs == B);
        mask_QN_AB = ((sum(mask_QN_A, 2) + sum(mask_QN_B, 2)) == 2);
        numCorners(i) = sum(mask_QN_AB);
    end
% ##### QuadPoints and FCNoints #####
elseif nargin == 4
    for i = 1:size(labels, 1)
        A = labels(i, 1);
        B = labels(i, 2);

        mask_TL_A = (TLs == A);
        mask_TL_B = (TLs == B);
        mask_TL_AB = ((sum(mask_TL_A, 2) + sum(mask_TL_B, 2)) == 2);
        numEdges(i) = sum(mask_TL_AB);

        mask_QN_A = (QNs == A);
        mask_QN_B = (QNs == B);
        mask_QN_AB = ((sum(mask_QN_A, 2) + sum(mask_QN_B, 2)) == 2);
        
        mask_FCN_A = (FCNs == A);
        mask_FCN_B = (FCNs == B);
        mask_FCN_AB = ((sum(mask_FCN_A, 2) + sum(mask_FCN_B, 2)) == 2);
        numCorners(i) = sum(mask_QN_AB) + sum(mask_FCN_AB)*2;
    end
end

% mask_TL = zeros(length(TLs),2);
% for i = 1:length(TLs)
%     mask_TL(i, :) = ismember([A, B], TLs(i,:));
% end
% numEdges = sum(sum(mask_TL,2) == 2)
% mask_QN = zeros(length(generalQNs),2);
% for i = 1:length(generalQNs)
%     mask_QN(i, :) = ismember([A, B], generalQNs(i,:));
% end
% numCorners = sum(sum(mask_QN,2) == 2)
    
% end



end