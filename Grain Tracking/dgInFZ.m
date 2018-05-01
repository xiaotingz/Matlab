function [RFvec, misA] = dgInFZ(g1, g2, O)
% ##########################################################################
% * Input: 
%     - g1, g2, size = [3, 3]
%         orientation matrixes
% * Output: 
%     - RFvex, size = [3,1]
%         the variant of the misorientation, as Rodrigues vector, that lies in the Fundamental Zone
%     - g_FZ, size = [3, 3]
%         the variant of the misorientation, as orientation matrix, that lies in the Fundamental Zone         
% ##########################################################################
    for j = 1:24
        for k = 1:24
            gg1 = O(:,:,j)*g1;
            gg2 = O(:,:,k)*g2;
            
            dg = gg1*gg2';
            misA = acosd(0.5*(trace(dg)-1));
% # # misA error check-----------------------------------------------------------------
% #             if ((0.5*(trace(dg)-1) <= 1) & (0.5*(trace(dg)-1) >= -1)):
% #                 misA = acosd(0.5*(trace(dg)-1))
% #             else:
% #                 print 0.5*(trace(dg)-1
% #                 print 'misA problem at i =', i
% #                 return 
% # # misA error check-----------------------------------------------------------------
            
            RFvec = [dg(2,3)-dg(3,2), dg(3,1)-dg(1,3), dg(1,2)-dg(2,1)] /(2*sind(misA)) * tand(misA/2);
            inFZ = ((all(RFvec >= 0)) & (RFvec(1) <= (sqrt(2) - 1)) & ...
                      (RFvec(1) >= RFvec(2)) & (RFvec(2) >= RFvec(3)) & (sum(RFvec) <= 1) & (sum(RFvec) > 0.0001) );
            if inFZ
                return 
            end
            dg = gg2*gg1';
            misA = acosd(0.5*(trace(dg)-1));
% # # misA error check-----------------------------------------------------------------
% #             if ((0.5*(trace(dg)-1) <= 1) & (0.5*(trace(dg)-1) >= -1)):
% #                 misA = acosd(0.5*(trace(dg)-1))
% #             else:
% #                 print 0.5*(trace(dg)-1
% #                 print 'misA problem at i =', i
% #                 return
% # # misA error check-----------------------------------------------------------------                
            
            RFvec = [dg(2,3)-dg(3,2), dg(3,1)-dg(1,3), dg(1,2)-dg(2,1)] /(2*sind(misA)) * tand(misA/2);
            inFZ = ((all(RFvec >= 0)) & (RFvec(1) <= (sqrt(2) - 1)) & ...
                      (RFvec(1) >= RFvec(2)) & (RFvec(2) >= RFvec(3)) & (sum(RFvec) <= 1) & (sum(RFvec) > 0.0001));
            if inFZ
                return 
            end
        end
    end
