function v4dist = spv4dist(b_1, b_2, bsym)
% ##########################################################################
% * Input
%     - b = [4, 4]
%           Grain boundary representation matrix as defined by Morawiec.
%           Following Yufeng & Morawiec's formula.
%     - bsym
%           Symmetry operators
% * see FIND_NEIGH.F90 and GBdist.ipynb
% ##########################################################################
v4dist = 5.0;
for ii = 1:24
    for jj  = 1:24
        ttt = bsym(:,:,ii) * b_2 * bsym(:,:,jj);
        tt = b_1 - ttt;
        compar = trace(tt * tt');
        if compar < v4dist
            v4dist = compar;
            bb = ttt;
            truei = ii;
            truej = jj;
            tran = 1;
        end
        tt = b_1 - ttt';
        compar = trace(tt * tt');
        if compar < v4dist
            v4dist = compar;
            bb = ttt';
            truei = ii;
            truej = jj;
            tran = 2;
        end
        
        ttt(4, 1:3) = - ttt(4, 1:3);
        ttt(1:3, 4) = - ttt(1:3, 4);
        tt = b_1 - ttt;
        compar = trace(tt * tt');
        if compar < v4dist
            v4dist = compar;
            bb = ttt;
            truei = ii;
            truej = jj;
            tran = 3;
        end
        tt = b_1 - ttt';
        compar = trace(tt * tt');
        if compar < v4dist
            v4dist = compar;
            bb = ttt';
            truei = ii;
            truej = jj;
            tran = 4;
        end
    end
end

end

