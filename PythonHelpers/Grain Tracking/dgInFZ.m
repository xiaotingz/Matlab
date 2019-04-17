function [rf_vec, misA, dg] = dgInFZ(g1, g2, O)
% ##########################################################################
% * Input: 
%     - g1, g2, size = [3, 3]
%         orientation matrixes
% * Output: 
%     - RFvec, size = [3,1]
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
            
            rf_vec = [dg(2,3)-dg(3,2), dg(3,1)-dg(1,3), dg(1,2)-dg(2,1)] /(2*sind(misA)) * tand(misA/2);
            in_fz = ((all(rf_vec >= 0)) & (rf_vec(1) <= (sqrt(2) - 1)) & ...
                      (rf_vec(1) >= rf_vec(2)) & (rf_vec(2) >= rf_vec(3)) & (sum(rf_vec) <= 1) & (sum(rf_vec) > 0.0001) );
            if in_fz
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
            
            rf_vec = [dg(2,3)-dg(3,2), dg(3,1)-dg(1,3), dg(1,2)-dg(2,1)] /(2*sind(misA)) * tand(misA/2);
            in_fz = ((all(rf_vec >= 0)) & (rf_vec(1) <= (sqrt(2) - 1)) & ...
                      (rf_vec(1) >= rf_vec(2)) & (rf_vec(2) >= rf_vec(3)) & (sum(rf_vec) <= 1) & (sum(rf_vec) > 0.0001));
            if in_fz
                return 
            end
        end
    end
